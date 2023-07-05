// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "domainpostoperator.hpp"

#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

DomainPostOperator::DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       mfem::ParFiniteElementSpace *nd_fespace,
                                       mfem::ParFiniteElementSpace *rt_fespace,
                                       int pa_order_threshold)
{
  if (nd_fespace)
  {
    // Construct ND mass matrix to compute the electric field energy integral as:
    //              E_elec = 1/2 Re{∫_Ω Dᴴ E dV} as (M_eps * e)ᴴ e.
    // Only the real part of the permeability contributes to the energy (imaginary part
    // cancels out in the inner product due to symmetry).
    constexpr int skip_zeros = 0;
    constexpr auto MatTypeEpsReal = MaterialPropertyType::PERMITTIVITY_REAL;
    constexpr auto MatTypeEpsImag = MaterialPropertyType::PERMITTIVITY_IMAG;
    MaterialPropertyCoefficient<MatTypeEpsReal> epsilon_func(mat_op);
    auto m_nd = std::make_unique<mfem::SymmetricBilinearForm>(nd_fespace);
    m_nd->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(epsilon_func));
    M_ND = utils::AssembleOperator(std::move(m_nd), true, pa_order_threshold, skip_zeros);
    D.SetSize(M_ND->Height());

    // Use the provided domain postprocessing indices to group for postprocessing bulk
    // dielectric loss.
    int attr_max = nd_fespace->GetParMesh()->attributes.Max();
    for (const auto &[idx, data] : iodata.domains.postpro.dielectric)
    {
      mfem::Array<int> attr_marker(attr_max);
      attr_marker = 0;
      for (auto attr : data.attributes)
      {
        attr_marker[attr - 1] = 1;
      }
      SumMatrixCoefficient epsilon_func_r(nd_fespace->GetParMesh()->SpaceDimension());
      SumMatrixCoefficient epsilon_func_i(nd_fespace->GetParMesh()->SpaceDimension());
      epsilon_func_r.AddCoefficient(
          std::make_unique<MaterialPropertyCoefficient<MatTypeEpsReal>>(mat_op),
          attr_marker);
      epsilon_func_i.AddCoefficient(
          std::make_unique<MaterialPropertyCoefficient<MatTypeEpsImag>>(mat_op, -1.0),
          attr_marker);
      auto mr_nd = std::make_unique<mfem::SymmetricBilinearForm>(nd_fespace);
      auto mi_nd = std::make_unique<mfem::SymmetricBilinearForm>(nd_fespace);
      mr_nd->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(epsilon_func_r));
      mi_nd->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(epsilon_func_i));
      M_NDi.emplace(
          idx, std::make_pair(utils::AssembleOperator(std::move(mr_nd), true,
                                                      pa_order_threshold, skip_zeros),
                              utils::AssembleOperator(std::move(mi_nd), true,
                                                      pa_order_threshold, skip_zeros)));
    }
  }

  if (rt_fespace)
  {
    // Construct RT mass matrix to compute the magnetic field energy integral as:
    //              E_mag = 1/2 Re{∫_Ω Bᴴ H dV} as (M_muinv * b)ᴴ b.
    constexpr int skip_zeros = 0;
    constexpr auto MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
    MaterialPropertyCoefficient<MatTypeMuInv> muinv_func(mat_op);
    auto m_rt = std::make_unique<mfem::SymmetricBilinearForm>(rt_fespace);
    m_rt->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(muinv_func));
    M_RT =
        utils::AssembleOperator(std::move(m_rt), true, pa_order_threshold - 1, skip_zeros);
    H.SetSize(M_RT->Height());
  }
}

double
DomainPostOperator::GetElectricFieldEnergy(const mfem::ParComplexGridFunction &E) const
{
  if (M_ND)
  {
    M_ND->Mult(E.real(), D);
    double res = mfem::InnerProduct(E.real(), D);
    M_ND->Mult(E.imag(), D);
    res += mfem::InnerProduct(E.imag(), D);
    Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for electric field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetElectricFieldEnergy(const mfem::ParGridFunction &E) const
{
  if (M_ND)
  {
    M_ND->Mult(E, D);
    double res = mfem::InnerProduct(E, D);
    Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for electric field energy calculation!");
  return 0.0;
}

double
DomainPostOperator::GetMagneticFieldEnergy(const mfem::ParComplexGridFunction &B) const
{
  if (M_RT)
  {
    M_RT->Mult(B.real(), H);
    double res = mfem::InnerProduct(B.real(), H);
    M_RT->Mult(B.imag(), H);
    res += mfem::InnerProduct(B.imag(), H);
    Mpi::GlobalSum(1, &res, B.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for magnetic field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetMagneticFieldEnergy(const mfem::ParGridFunction &B) const
{
  if (M_RT)
  {
    M_RT->Mult(B, H);
    double res = mfem::InnerProduct(B, H);
    Mpi::GlobalSum(1, &res, B.ParFESpace()->GetComm());
    return 0.5 * res;
  }
  MFEM_ABORT(
      "Domain postprocessing is not configured for magnetic field energy calculation!");
  return 0.0;
}

double DomainPostOperator::GetDomainElectricFieldEnergy(
    int idx, const mfem::ParComplexGridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = M_NDi.find(idx);
  MFEM_VERIFY(it != M_NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  it->second.first->Mult(E.real(), D);
  double res = mfem::InnerProduct(E.real(), D);
  it->second.first->Mult(E.imag(), D);
  res += mfem::InnerProduct(E.imag(), D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double
DomainPostOperator::GetDomainElectricFieldEnergy(int idx,
                                                 const mfem::ParGridFunction &E) const
{
  auto it = M_NDi.find(idx);
  MFEM_VERIFY(it != M_NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  it->second.first->Mult(E, D);
  double res = mfem::InnerProduct(E, D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double DomainPostOperator::GetDomainElectricFieldEnergyLoss(
    int idx, const mfem::ParComplexGridFunction &E) const
{
  // Compute the electric field energy integral for only a portion of the domain.
  auto it = M_NDi.find(idx);
  MFEM_VERIFY(it != M_NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  it->second.second->Mult(E.real(), D);
  double res = mfem::InnerProduct(E.real(), D);
  it->second.second->Mult(E.imag(), D);
  res += mfem::InnerProduct(E.imag(), D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

double
DomainPostOperator::GetDomainElectricFieldEnergyLoss(int idx,
                                                     const mfem::ParGridFunction &E) const
{
  auto it = M_NDi.find(idx);
  MFEM_VERIFY(it != M_NDi.end(),
              "Invalid domain index when postprocessing bulk dielectric loss!");
  it->second.second->Mult(E, D);
  double res = mfem::InnerProduct(E, D);
  Mpi::GlobalSum(1, &res, E.ParFESpace()->GetComm());
  return 0.5 * res;
}

}  // namespace palace
