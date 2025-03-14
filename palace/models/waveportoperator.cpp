// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "waveportoperator.hpp"

#include <array>
#include <unordered_map>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/arpack.hpp"
#include "linalg/iterative.hpp"
#include "linalg/mumps.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

void GetEssentialTrueDofs(mfem::ParGridFunction &E0t, mfem::ParGridFunction &E0n,
                          mfem::ParGridFunction &port_E0t, mfem::ParGridFunction &port_E0n,
                          mfem::ParTransferMap &port_nd_transfer,
                          mfem::ParTransferMap &port_h1_transfer,
                          const mfem::Array<int> &dbc_attr,
                          mfem::Array<int> &port_nd_dbc_tdof_list,
                          mfem::Array<int> &port_h1_dbc_tdof_list)
{
  auto &nd_fespace = *E0t.ParFESpace();
  auto &h1_fespace = *E0n.ParFESpace();
  auto &port_nd_fespace = *port_E0t.ParFESpace();
  auto &port_h1_fespace = *port_E0n.ParFESpace();
  const auto &mesh = *nd_fespace.GetParMesh();

  mfem::Array<int> dbc_marker, nd_dbc_tdof_list, h1_dbc_tdof_list;
  mesh::AttrToMarker(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0, dbc_attr,
                     dbc_marker);
  nd_fespace.GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);

  Vector tE0t(nd_fespace.GetTrueVSize()), tE0n(h1_fespace.GetTrueVSize());
  tE0t = 0.0;
  tE0n = 0.0;
  linalg::SetSubVector(tE0t, nd_dbc_tdof_list, 1.0);
  linalg::SetSubVector(tE0n, h1_dbc_tdof_list, 1.0);
  E0t.SetFromTrueDofs(tE0t);
  E0n.SetFromTrueDofs(tE0n);
  port_nd_transfer.Transfer(E0t, port_E0t);
  port_h1_transfer.Transfer(E0n, port_E0n);

  Vector port_tE0t(port_nd_fespace.GetTrueVSize()),
      port_tE0n(port_h1_fespace.GetTrueVSize());
  port_E0t.ParallelProject(port_tE0t);
  port_E0n.ParallelProject(port_tE0n);
  for (int i = 0; i < port_tE0t.Size(); i++)
  {
    if (port_tE0t[i] != 0.0)
    {
      port_nd_dbc_tdof_list.Append(i);
    }
  }
  for (int i = 0; i < port_tE0n.Size(); i++)
  {
    if (port_tE0n[i] != 0.0)
    {
      port_h1_dbc_tdof_list.Append(i);
    }
  }
}

void GetInitialSpace(const mfem::ParFiniteElementSpace &nd_fespace,
                     const mfem::ParFiniteElementSpace &h1_fespace,
                     const mfem::Array<int> &nd_dbc_tdof_list,
                     const mfem::Array<int> &h1_dbc_tdof_list, ComplexVector &v)
{
  // Initial space chosen as such that B v₀ = y₀, with y₀ = [y₀ₜ, 0, ... 0]ᵀ ⟂ null(A)
  // (with Aₜₜ nonsingular). See Lee, Sun, and Cendes, 1991 for reference.
  // Note: When the eigenvalue solver uses a standard ℓ²-inner product instead of B-inner
  // product (since we use a general non-Hermitian solver due to complex symmetric B), then
  // we just use v0 = y0 directly.
  const int nd_size = nd_fespace.GetTrueVSize(), h1_size = h1_fespace.GetTrueVSize();
  v.SetSize(nd_size + h1_size);
  // linalg::SetRandomReal(nd_fespace.GetComm(), v);
  v = std::complex<double>(1.0, 0.0);
  linalg::SetSubVector(v, nd_dbc_tdof_list, 0.0);
  linalg::SetSubVector(v, nd_size, nd_size + h1_size, 0.0);
}

constexpr bool skip_zeros = false;

std::unique_ptr<mfem::HypreParMatrix> GetBtt(const MaterialOperator &mat_op,
                                             const FiniteElementSpace &nd_fespace)
{
  // Mass matrix: Bₜₜ = (μ⁻¹ u, v).
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  return ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
}

std::unique_ptr<mfem::HypreParMatrix> GetBtn(const MaterialOperator &mat_op,
                                             const FiniteElementSpace &nd_fespace,
                                             const FiniteElementSpace &h1_fespace)
{
  // Mass matrix: Bₜₙ = (μ⁻¹ ∇ₜ u, v).
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  BilinearForm btn(h1_fespace, nd_fespace);
  btn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  return ParOperator(btn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
      .StealParallelAssemble();
}

std::array<std::unique_ptr<mfem::HypreParMatrix>, 2>
GetBnn(const MaterialOperator &mat_op, const FiniteElementSpace &h1_fespace,
       const mfem::Vector &normal, double omega)
{
  // Mass matrix: Bₙₙ = (μ⁻¹ ∇ₜ u, ∇ₜ v) - ω² (ε u, v).
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  MaterialPropertyCoefficient epsilon_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetPermittivityReal(), -omega * omega);
  epsilon_func.NormalProjectedCoefficient(normal);
  BilinearForm bnnr(h1_fespace);
  bnnr.AddDomainIntegrator<DiffusionMassIntegrator>(muinv_func, epsilon_func);

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {ParOperator(bnnr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble(),
            nullptr};
  }
  MaterialPropertyCoefficient negepstandelta_func(
      mat_op.GetBdrAttributeToMaterial(), mat_op.GetPermittivityImag(), -omega * omega);
  negepstandelta_func.NormalProjectedCoefficient(normal);
  BilinearForm bnni(h1_fespace);
  bnni.AddDomainIntegrator<MassIntegrator>(negepstandelta_func);
  return {ParOperator(bnnr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble(),
          ParOperator(bnni.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble()};
}

std::array<std::unique_ptr<mfem::HypreParMatrix>, 2>
GetAtt(const MaterialOperator &mat_op, const FiniteElementSpace &nd_fespace,
       const mfem::Vector &normal, double omega, double theta2)
{
  // Stiffness matrix (shifted):
  //     Aₜₜ = 1/Θ² [(μ⁻¹ ∇ₜ x u, ∇ₜ x v) - ω² (ε u, v)] + (μ⁻¹ u, v).
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability(), 1.0 / theta2);
  muinv_func.NormalProjectedCoefficient(normal);
  MaterialPropertyCoefficient epsilon_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetPermittivityReal(),
                                           -omega * omega / theta2);
  epsilon_func.AddCoefficient(mat_op.GetBdrAttributeToMaterial(),
                              mat_op.GetInvPermeability());
  BilinearForm attr(nd_fespace);
  attr.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_func, epsilon_func);

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {ParOperator(attr.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
            nullptr};
  }
  MaterialPropertyCoefficient negepstandelta_func(mat_op.GetBdrAttributeToMaterial(),
                                                  mat_op.GetPermittivityImag(),
                                                  -omega * omega / theta2);
  BilinearForm atti(nd_fespace);
  atti.AddDomainIntegrator<VectorFEMassIntegrator>(negepstandelta_func);
  return {ParOperator(attr.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          ParOperator(atti.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble()};
}

std::array<std::unique_ptr<mfem::HypreParMatrix>, 4>
GetSystemMatrices(mfem::HypreParMatrix *Btt, mfem::HypreParMatrix *Btn,
                  mfem::HypreParMatrix *BtnT, mfem::HypreParMatrix *Bnnr,
                  mfem::HypreParMatrix *Bnni, mfem::HypreParMatrix *Attr,
                  mfem::HypreParMatrix *Atti, mfem::HypreParMatrix *Dtt,
                  const mfem::Array<int> &dbc_tdof_list)
{
  // Construct the 2x2 block matrices for the eigenvalue problem K e = λ M e.
  mfem::Array2D<mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Btt;
  blocks(0, 1) = Btn;
  blocks(1, 0) = BtnT;
  blocks(1, 1) = Bnnr;
  std::unique_ptr<mfem::HypreParMatrix> Kr(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Ki;
  if (Bnni)
  {
    blocks = nullptr;
    blocks(0, 0) = Dtt;
    blocks(1, 1) = Bnni;
    Ki.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  blocks(0, 0) = Attr;  // Att is already shifted (= Bₜₜ + 1/Θ² Aₜₜ)
  blocks(0, 1) = Btn;
  blocks(1, 0) = BtnT;
  blocks(1, 1) = Bnnr;
  std::unique_ptr<mfem::HypreParMatrix> Mr(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Mi;
  if (Atti)
  {
    MFEM_VERIFY(Bnni, "Both imaginary parts should exist for wave port eigenvalue solver!");
    blocks = nullptr;
    blocks(0, 0) = Atti;  // Att is already shifted (= Bₜₜ + 1/Θ² Aₜₜ)
    blocks(1, 1) = Bnni;
    Mi.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs not associated with this wave port or constrained by
  // Dirichlet BCs.
  Kr->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  if (Ki)
  {
    Ki->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }
  Mr->EliminateBC(dbc_tdof_list, Operator::DIAG_ONE);
  if (Mi)
  {
    Mi->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(Kr), std::move(Ki), std::move(Mr), std::move(Mi)};
}

void NormalizeWithSign(const mfem::ParGridFunction &S0t, mfem::ParComplexGridFunction &E0t,
                       mfem::ParComplexGridFunction &E0n, mfem::LinearForm &sr,
                       mfem::LinearForm &si)
{
  // Normalize grid functions to a chosen polarization direction and unit power, |E x H⋆| ⋅
  // n, integrated over the port surface (+n is the direction of propagation). The n x H
  // coefficients are updated implicitly as the only store references to the Et, En grid
  // functions. We choose a (rather arbitrary) sign constraint to at least make results for
  // the same port consistent between frequencies/meshes.

  // |E x H⋆| ⋅ n = |E ⋅ (-n x H⋆)|
  double sign = sr * S0t;
  std::complex<double> dot(-(sr * E0t.real()) - (si * E0t.imag()),
                           -(sr * E0t.imag()) + (si * E0t.real()));
  std::array<double, 3> data = {sign, dot.real(), dot.imag()};
  Mpi::GlobalSum(3, data.data(), S0t.ParFESpace()->GetComm());
  sign = (data[0] < 0.0) ? -1.0 : 1.0;
  dot = {data[1], data[2]};

  double scale = sign / std::sqrt(std::abs(dot));
  E0t.real() *= scale;  // Updates the n x H coefficients depending on Et, En too
  E0t.imag() *= scale;
  E0n.real() *= scale;
  E0n.imag() *= scale;
  sr *= scale;  // Update linear forms for postprocessing
  si *= scale;

  // This parallel communication is not required since wave port boundaries are true
  // one-sided boundaries.
  // port_E0t->real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
  // port_E0t->imag().ExchangeFaceNbrData();  // for n x H coefficients evaluation
  // port_E0n->real().ExchangeFaceNbrData();
  // port_E0n->imag().ExchangeFaceNbrData();
}

// Helper for BdrSubmeshEVectorCoefficient and BdrSubmeshHVectorCoefficient.
enum class ValueType
{
  REAL,
  IMAG
};

// Return as a vector coefficient the boundary mode electric field.
template <ValueType Type>
class BdrSubmeshEVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const mfem::ParComplexGridFunction &Et, &En;
  const mfem::ParSubMesh &submesh;
  const std::unordered_map<int, int> &submesh_parent_elems;
  mfem::IsoparametricTransformation T_loc;

public:
  BdrSubmeshEVectorCoefficient(const mfem::ParComplexGridFunction &Et,
                               const mfem::ParComplexGridFunction &En,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems)
    : mfem::VectorCoefficient(Et.real().VectorDim()), Et(Et), En(En), submesh(submesh),
      submesh_parent_elems(submesh_parent_elems)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Always do the GridFunction evaluation in the submesh.
    mfem::ElementTransformation *T_submesh = nullptr;
    if (T.mesh == submesh.GetParent())
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                  "BdrSubmeshEVectorCoefficient requires ElementType::BDR_ELEMENT when not "
                  "used on a SubMesh!");
      auto it = submesh_parent_elems.find(T.ElementNo);
      if (it == submesh_parent_elems.end())
      {
        // Just return zero for a parent boundary element not in the submesh.
        V.SetSize(vdim);
        V = 0.0;
        return;
      }
      else
      {
        submesh.GetElementTransformation(it->second, &T_loc);
        T_loc.SetIntPoint(&ip);
        T_submesh = &T_loc;
      }
    }
    else if (T.mesh == &submesh)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                  "BdrSubmeshEVectorCoefficient requires ElementType::ELEMENT when used on "
                  "a SubMesh!");
      T_submesh = &T;
    }
    else
    {
      MFEM_ABORT("Invalid mesh for BdrSubmeshEVectorCoefficient!");
    }

    // Compute Eₜ + n ⋅ Eₙ .
    mfem::Vector U, nor;
    BdrGridFunctionCoefficient::GetNormal(*T_submesh, nor);
    if constexpr (Type == ValueType::REAL)
    {
      Et.real().GetVectorValue(*T_submesh, ip, V);
      auto Vn = En.real().GetValue(*T_submesh, ip);
      V.Add(Vn, nor);
    }
    else
    {
      Et.imag().GetVectorValue(*T_submesh, ip, V);
      auto Vn = En.imag().GetValue(*T_submesh, ip);
      V.Add(Vn, nor);
    }
  }
};

// Computes boundary mode n x H, where +n is the direction of wave propagation: n x H =
// -1/(iωμ) (ikₙ Eₜ + ∇ₜ Eₙ), using the tangential and normal electric field component grid
// functions evaluated on the (single-sided) boundary element.
template <ValueType Type>
class BdrSubmeshHVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const mfem::ParComplexGridFunction &Et, &En;
  const MaterialOperator &mat_op;
  const mfem::ParSubMesh &submesh;
  const std::unordered_map<int, int> &submesh_parent_elems;
  mfem::IsoparametricTransformation T_loc;
  std::complex<double> kn;
  double omega;

public:
  BdrSubmeshHVectorCoefficient(const mfem::ParComplexGridFunction &Et,
                               const mfem::ParComplexGridFunction &En,
                               const MaterialOperator &mat_op,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems,
                               std::complex<double> kn, double omega)
    : mfem::VectorCoefficient(Et.real().VectorDim()), Et(Et), En(En), mat_op(mat_op),
      submesh(submesh), submesh_parent_elems(submesh_parent_elems), kn(kn), omega(omega)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Always do the GridFunction evaluation in the submesh.
    mfem::ElementTransformation *T_submesh = nullptr;
    if (T.mesh == submesh.GetParent())
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                  "BdrSubmeshHVectorCoefficient requires ElementType::BDR_ELEMENT when not "
                  "used on a SubMesh!");
      auto it = submesh_parent_elems.find(T.ElementNo);
      if (it == submesh_parent_elems.end())
      {
        // Just return zero for a parent boundary element not in the submesh.
        V.SetSize(vdim);
        V = 0.0;
        return;
      }
      else
      {
        submesh.GetElementTransformation(it->second, &T_loc);
        T_loc.SetIntPoint(&ip);
        T_submesh = &T_loc;
      }
    }
    else if (T.mesh == &submesh)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                  "BdrSubmeshHVectorCoefficient requires ElementType::ELEMENT when used on "
                  "a SubMesh!");
      T_submesh = &T;
    }
    else
    {
      MFEM_ABORT("Invalid mesh for BdrSubmeshHVectorCoefficient!");
    }

    // Get the attribute in the neighboring domain element of the parent mesh.
    int attr = [&T, this]()
    {
      int i = -1, o, iel1, iel2;
      if (T.mesh == submesh.GetParent())
      {
        MFEM_ASSERT(
            T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
            "BdrSubmeshHVectorCoefficient requires ElementType::BDR_ELEMENT when not "
            "used on a SubMesh!");
        T.mesh->GetBdrElementFace(T.ElementNo, &i, &o);
      }
      else if (T.mesh == &submesh)
      {
        MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                    "BdrSubmeshHVectorCoefficient requires ElementType::ELEMENT when used "
                    "on a SubMesh!");
        submesh.GetParent()->GetBdrElementFace(submesh.GetParentElementIDMap()[T.ElementNo],
                                               &i, &o);
      }
      else
      {
        MFEM_ABORT("Invalid mesh for BdrSubmeshHVectorCoefficient!");
      }
      submesh.GetParent()->GetFaceElements(i, &iel1, &iel2);
      return submesh.GetParent()->GetAttribute(iel1);
    }();

    // Compute Re/Im{-1/i (ikₙ Eₜ + ∇ₜ Eₙ)} (t-gradient evaluated in boundary element).
    mfem::Vector U;
    if constexpr (Type == ValueType::REAL)
    {
      Et.real().GetVectorValue(*T_submesh, ip, U);
      U *= -kn.real();

      mfem::Vector dU;
      En.imag().GetGradient(*T_submesh, dU);
      U -= dU;
    }
    else
    {
      Et.imag().GetVectorValue(*T_submesh, ip, U);
      U *= -kn.real();

      mfem::Vector dU;
      En.real().GetGradient(*T_submesh, dU);
      U += dU;
    }

    // Scale by 1/(ωμ) with μ evaluated in the neighboring element.
    V.SetSize(U.Size());
    mat_op.GetInvPermeability(attr).Mult(U, V);
    V *= (1.0 / omega);
  }
};

}  // namespace

WavePortData::WavePortData(const config::WavePortData &data,
                           const config::SolverData &solver, const MaterialOperator &mat_op,
                           mfem::ParFiniteElementSpace &nd_fespace,
                           mfem::ParFiniteElementSpace &h1_fespace,
                           const mfem::Array<int> &dbc_attr)
  : mat_op(mat_op)
{
  mode_idx = data.mode_idx;
  d_offset = data.d_offset;
  excitation = data.excitation;
  kn0 = 0.0;
  omega0 = 0.0;

  // Construct the SubMesh.
  MFEM_VERIFY(!data.attributes.empty(), "Wave port boundary found with no attributes!");
  const auto &mesh = *nd_fespace.GetParMesh();
  attr_list.Append(data.attributes.data(), data.attributes.size());
  port_mesh = std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromBoundary(mesh, attr_list)));
  port_normal = mesh::GetSurfaceNormal(*port_mesh);

  port_nd_fec = std::make_unique<mfem::ND_FECollection>(nd_fespace.GetMaxElementOrder(),
                                                        port_mesh->Dimension());
  port_h1_fec = std::make_unique<mfem::H1_FECollection>(h1_fespace.GetMaxElementOrder(),
                                                        port_mesh->Dimension());
  port_nd_fespace = std::make_unique<FiniteElementSpace>(*port_mesh, port_nd_fec.get());
  port_h1_fespace = std::make_unique<FiniteElementSpace>(*port_mesh, port_h1_fec.get());

  mfem::ParGridFunction E0t(&nd_fespace);
  mfem::ParGridFunction E0n(&h1_fespace);
  port_E0t = std::make_unique<mfem::ParComplexGridFunction>(&port_nd_fespace->Get());
  port_E0n = std::make_unique<mfem::ParComplexGridFunction>(&port_h1_fespace->Get());

  port_nd_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0t, port_E0t->real()));
  port_h1_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0n, port_E0n->real()));

  // Construct mapping from parent (boundary) element indices to submesh (domain)
  // elements.
  {
    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    const mfem::Array<int> &parent_elems = port_submesh.GetParentElementIDMap();
    for (int i = 0; i < parent_elems.Size(); i++)
    {
      submesh_parent_elems[parent_elems[i]] = i;
    }
  }

  // The operators for the generalized eigenvalue problem are:
  //                [Aₜₜ  0] [eₜ]  = -kₙ² [Bₜₜ   Bₜₙ] [eₜ]
  //                [0   0] [eₙ]        [Bₜₙᵀ  Bₙₙ] [eₙ]
  // for the wave port of the given index. The transformed variables are related to the true
  // field by Eₜ = eₜ/kₙ and Eₙ = ieₙ. This is solved on the global mesh so the result is a
  // grid function over the entire space, not just the port boundary (so that it can be
  // queried from functions which use the global mesh).
  //
  // We will actually solve the shifted problem K e = λ M e, or:
  //                [Bₜₜ   Bₜₙ] [eₜ]   =  λ [Bₜₜ + 1/Θ² Aₜₜ  Bₜₙ] [eₜ]
  //                [Bₜₙᵀ  Bₙₙ] [eₙ]       [Bₜₙᵀ          Bₙₙ] [eₙ] .
  // Here we have λ = Θ²/(Θ²-kₙ²), where Θ² bounds the maximum kₙ² and is taken as Θ² =
  // ω² μₘₐₓ εₘₐₓ over the entire simulation domain.
  // Reference: Lee, Sun, and Cendes, Full-wave analysis of dielectric waveguides using
  //            tangential vector finite elements, IEEE Trans. Microwave Theory Tech.
  //            (1991).
  double c_min = mat_op.GetLightSpeedMin().Min();
  MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
              "Invalid material speed of light detected in WavePortOperator!");
  mu_eps_max = 1.0 / (c_min * c_min);
  Btt = GetBtt(mat_op, *port_nd_fespace);
  Btn = GetBtn(mat_op, *port_nd_fespace, *port_h1_fespace);
  BtnT.reset(Btn->Transpose());
  {
    // The HyperParMatrix constructor from a SparseMatrix on each process does not copy the
    // SparseMatrix data, so we deep copy the entire data structure manually.
    mfem::Vector d(Btt->Height());
    d.UseDevice(false);  // SparseMatrix constructor uses Vector on host
    d = 0.0;
    mfem::SparseMatrix diag(d);
    Dtt = std::make_unique<mfem::HypreParMatrix>(mfem::HypreParMatrix(
        Btt->GetComm(), Btt->GetGlobalNumRows(), Btt->GetRowStarts(), &diag));
  }

  // Extract Dirichlet BC true dofs for the port FE spaces.
  mfem::Array<int> port_nd_dbc_tdof_list, port_h1_dbc_tdof_list;
  {
    GetEssentialTrueDofs(E0t, E0n, port_E0t->real(), port_E0n->real(), *port_nd_transfer,
                         *port_h1_transfer, dbc_attr, port_nd_dbc_tdof_list,
                         port_h1_dbc_tdof_list);
    int nd_tdof_offset = port_nd_fespace->GetTrueVSize();
    port_dbc_tdof_list.Reserve(port_nd_dbc_tdof_list.Size() + port_h1_dbc_tdof_list.Size());
    for (auto tdof : port_nd_dbc_tdof_list)
    {
      port_dbc_tdof_list.Append(tdof);
    }
    for (auto tdof : port_h1_dbc_tdof_list)
    {
      port_dbc_tdof_list.Append(tdof + nd_tdof_offset);
    }
  }

  // Create vector for initial space for eigenvalue solves (for nullspace of [Aₜₜ  0]
  //                                                                         [0   0]).
  GetInitialSpace(*port_nd_fespace, *port_h1_fespace, port_nd_dbc_tdof_list,
                  port_h1_dbc_tdof_list, v0);
  e0.SetSize(v0.Size());
  e0t.SetSize(port_nd_fespace->GetTrueVSize());
  e0n.SetSize(port_h1_fespace->GetTrueVSize());

  // Configure a communicator for the processes which have elements for this port.
  MPI_Comm comm = nd_fespace.GetComm();
  int color = (port_nd_fespace->GetVSize() > 0 || port_h1_fespace->GetVSize() > 0)
                  ? 0
                  : MPI_UNDEFINED;
  MPI_Comm_split(comm, color, Mpi::Rank(comm), &port_comm);
  MFEM_VERIFY((color == 0 && port_comm != MPI_COMM_NULL) ||
                  (color == MPI_UNDEFINED && port_comm == MPI_COMM_NULL),
              "Unexpected error splitting communicator for wave port boundaries!");
  port_root = (color == MPI_UNDEFINED) ? Mpi::Size(comm) : Mpi::Rank(comm);
  Mpi::GlobalMin(1, &port_root, comm);
  MFEM_VERIFY(port_root < Mpi::Size(comm), "No root process found for port!");

  // Configure the eigenvalue problem solver. As for the full 3D case, the system matrices
  // are in general complex and symmetric. We supply the operators to the solver in
  // shift-inverted form and handle the back-transformation externally.
  if (port_comm != MPI_COMM_NULL)
  {
    // Define the linear solver to be used for solving systems associated with the
    // generalized eigenvalue problem.
    constexpr int ksp_print = 0;
    constexpr double ksp_tol = 1.0e-8;
    constexpr double ksp_max_it = 30;
    auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(port_comm, ksp_print);
    gmres->SetInitialGuess(false);
    gmres->SetRelTol(ksp_tol);
    gmres->SetMaxIter(ksp_max_it);
    gmres->SetRestartDim(ksp_max_it);
    // gmres->SetPrecSide(GmresSolver<ComplexOperator>::PrecSide::RIGHT);

    config::LinearSolverData::Type pc_type = solver.linear.type;
    if (pc_type == config::LinearSolverData::Type::SUPERLU)
    {
#if !defined(MFEM_USE_SUPERLU)
      MFEM_ABORT("Solver was not built with SuperLU_DIST support, please choose a "
                 "different solver!");
#endif
    }
    else if (pc_type == config::LinearSolverData::Type::STRUMPACK ||
             pc_type == config::LinearSolverData::Type::STRUMPACK_MP)
    {
#if !defined(MFEM_USE_STRUMPACK)
      MFEM_ABORT("Solver was not built with STRUMPACK support, please choose a "
                 "different solver!");
#endif
    }
    else if (pc_type == config::LinearSolverData::Type::MUMPS)
    {
#if !defined(MFEM_USE_MUMPS)
      MFEM_ABORT("Solver was not built with MUMPS support, please choose a "
                 "different solver!");
#endif
    }
    else  // Default choice
    {
#if defined(MFEM_USE_SUPERLU)
      pc_type = config::LinearSolverData::Type::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
      pc_type = config::LinearSolverData::Type::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
      pc_type = config::LinearSolverData::Type::MUMPS;
#else
#error "Wave port solver requires building with SuperLU_DIST, STRUMPACK, or MUMPS!"
#endif
    }
    auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
        [&]() -> std::unique_ptr<mfem::Solver>
        {
          if (pc_type == config::LinearSolverData::Type::SUPERLU)
          {
#if defined(MFEM_USE_SUPERLU)
            auto slu = std::make_unique<SuperLUSolver>(
                port_comm, config::LinearSolverData::SymFactType::DEFAULT, false,
                ksp_print - 1);
            slu->GetSolver().SetColumnPermutation(mfem::superlu::MMD_AT_PLUS_A);
            return slu;
#endif
          }
          else if (pc_type == config::LinearSolverData::Type::STRUMPACK)
          {
#if defined(MFEM_USE_STRUMPACK)
            auto strumpack = std::make_unique<StrumpackSolver>(
                port_comm, config::LinearSolverData::SymFactType::DEFAULT,
                config::LinearSolverData::CompressionType::NONE, 0.0, 0, 0, ksp_print - 1);
            strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::AMD);
            return strumpack;
#endif
          }
          else if (pc_type == config::LinearSolverData::Type::MUMPS)
          {
#if defined(MFEM_USE_MUMPS)
            auto mumps = std::make_unique<MumpsSolver>(
                port_comm, mfem::MUMPSSolver::SYMMETRIC_INDEFINITE,
                config::LinearSolverData::SymFactType::DEFAULT, 0.0, ksp_print - 1);
            mumps->SetReorderingStrategy(mfem::MUMPSSolver::AMD);
            return mumps;
#endif
          }
          return {};
        }());
    pc->SetSaveAssembled(false);
    ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));

    // Define the eigenvalue solver.
    constexpr int print = 0;
    config::EigenSolverData::Type type = solver.eigenmode.type;
    if (type == config::EigenSolverData::Type::SLEPC)
    {
#if !defined(PALACE_WITH_SLEPC)
      MFEM_ABORT("Solver was not built with SLEPc support, please choose a "
                 "different solver!");
#endif
    }
    else if (type == config::EigenSolverData::Type::ARPACK)
    {
#if !defined(PALACE_WITH_ARPACK)
      MFEM_ABORT("Solver was not built with ARPACK support, please choose a "
                 "different solver!");
#endif
    }
    else  // Default choice
    {
#if defined(PALACE_WITH_SLEPC)
      type = config::EigenSolverData::Type::SLEPC;
#elif defined(PALACE_WITH_ARPACK)
      type = config::EigenSolverData::Type::ARPACK;
#else
#error "Wave port solver requires building with ARPACK or SLEPc!"
#endif
    }
    if (type == config::EigenSolverData::Type::ARPACK)
    {
#if defined(PALACE_WITH_ARPACK)
      eigen = std::make_unique<arpack::ArpackEPSSolver>(port_comm, print);
#endif
    }
    else  // config::EigenSolverData::Type::SLEPC
    {
#if defined(PALACE_WITH_SLEPC)
      auto slepc = std::make_unique<slepc::SlepcEPSSolver>(port_comm, print);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
      eigen = std::move(slepc);
#endif
    }
    constexpr double tol = 1.0e-6;
    eigen->SetNumModes(mode_idx, std::max(2 * mode_idx + 1, 5));
    eigen->SetTol(tol);
    eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_MAGNITUDE);
    eigen->SetLinearSolver(*ksp);
  }

  // Configure port mode sign convention: 1ᵀ Re{-n x H} >= 0 on the "upper-right quadrant"
  // of the wave port boundary, in order to deal with symmetry effectively.
  {
    Vector bbmin, bbmax;
    port_mesh->Get().GetBoundingBox(bbmin, bbmax);
    const int dim = port_mesh->SpaceDimension();

    double la = 0.0, lb = 0.0;
    int da = -1, db = -1;
    for (int d = 0; d < dim; d++)
    {
      double diff = bbmax(d) - bbmin(d);
      if (diff > la)
      {
        lb = la;
        la = diff;
        db = da;
        da = d;
      }
      else if (diff > lb)
      {
        lb = diff;
        db = d;
      }
    }
    MFEM_VERIFY(da >= 0 && db >= 0 && da != db,
                "Unexpected wave port geometry for normalization!");
    double ca = 0.5 * (bbmax[da] + bbmin[da]), cb = 0.5 * (bbmax[db] + bbmin[db]);

    auto TDirection = [da, db, ca, cb, dim](const Vector &x, Vector &f)
    {
      MFEM_ASSERT(x.Size() == dim,
                  "Invalid dimension mismatch for wave port mode normalization!");
      f.SetSize(dim);
      if (x[da] >= ca && x[db] >= cb)
      {
        f = 1.0;
      }
      else
      {
        f = 0.0;
      }
    };
    mfem::VectorFunctionCoefficient tfunc(dim, TDirection);
    port_S0t = std::make_unique<mfem::ParGridFunction>(&port_nd_fespace->Get());
    port_S0t->ProjectCoefficient(tfunc);
  }
}

WavePortData::~WavePortData()
{
  // Free the solvers before the communicator on which they are based.
  ksp.reset();
  eigen.reset();
  if (port_comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&port_comm);
  }
}

void WavePortData::Initialize(double omega)
{
  if (omega == omega0)
  {
    return;
  }

  // Construct matrices and solve the generalized eigenvalue problem for the desired wave
  // port mode.
  double theta2 = mu_eps_max * omega * omega;
  std::unique_ptr<ComplexOperator> K, M;
  {
    auto [Bnnr, Bnni] = GetBnn(mat_op, *port_h1_fespace, port_normal, omega);
    auto [Attr, Atti] = GetAtt(mat_op, *port_nd_fespace, port_normal, omega, theta2);
    auto [Kr, Ki, Mr, Mi] =
        GetSystemMatrices(Btt.get(), Btn.get(), BtnT.get(), Bnnr.get(), Bnni.get(),
                          Attr.get(), Atti.get(), Dtt.get(), port_dbc_tdof_list);
    K = std::make_unique<ComplexWrapperOperator>(std::move(Kr), std::move(Ki));
    M = std::make_unique<ComplexWrapperOperator>(std::move(Mr), std::move(Mi));
  }

  // Configure and solve the eigenvalue problem for the desired boundary mode.
  std::complex<double> lambda;
  if (port_comm != MPI_COMM_NULL)
  {
    ComplexWrapperOperator P(M->Real(), nullptr);  // Non-owning constructor
    ksp->SetOperators(*M, P);
    eigen->SetOperators(*K, *M, EigenvalueSolver::ScaleType::NONE);
    eigen->SetInitialSpace(v0);
    int num_conv = eigen->Solve();
    MFEM_VERIFY(num_conv >= mode_idx, "Wave port eigensolver did not converge!");
    lambda = eigen->GetEigenvalue(mode_idx - 1);
    // Mpi::Print(port_comm, " ... Wave port eigensolver error = {} (bkwd), {} (abs)\n",
    //            eigen->GetError(mode_idx - 1, EigenvalueSolver::ErrorType::BACKWARD),
    //            eigen->GetError(mode_idx - 1, EigenvalueSolver::ErrorType::ABSOLUTE));
  }
  Mpi::Broadcast(1, &lambda, port_root, port_mesh->GetComm());

  // Extract the eigenmode solution and postprocess. The extracted eigenvalue is λ =
  // Θ² / (Θ² - kₙ²).
  MFEM_VERIFY(lambda.real() > 1.0 / (1.0 - 1.0e-2),
              "Computed wave port mode is or is very close to being evanescent "
                  << "(λ = " << lambda << ")!");
  kn0 = std::sqrt(theta2 - theta2 / lambda);
  omega0 = omega;

  // Separate the computed field out into eₜ and eₙ and and transform back to true
  // electric field variables: Eₜ = eₜ/kₙ and Eₙ = ieₙ.
  if (port_comm != MPI_COMM_NULL)
  {
    Vector e0tr, e0ti, e0nr, e0ni;
    eigen->GetEigenvector(mode_idx - 1, e0);
    e0tr.MakeRef(e0.Real(), 0, e0t.Size());
    e0nr.MakeRef(e0.Real(), e0t.Size(), e0n.Size());
    e0ti.MakeRef(e0.Imag(), 0, e0t.Size());
    e0ni.MakeRef(e0.Imag(), e0t.Size(), e0n.Size());
    e0t.Real() = e0tr;
    e0t.Imag() = e0ti;
    e0n.Real() = e0nr;
    e0n.Imag() = e0ni;
    e0t *= 1.0 / kn0;
    e0n *= 1i;
  }
  else
  {
    MFEM_ASSERT(e0.Size() == 0 && e0t.Size() == 0 && e0n.Size() == 0,
                "Unexpected non-empty port FE space in wave port boundary mode solve!");
  }
  port_E0t->real().SetFromTrueDofs(e0t.Real());  // Parallel distribute
  port_E0t->imag().SetFromTrueDofs(e0t.Imag());
  port_E0n->real().SetFromTrueDofs(e0n.Real());
  port_E0n->imag().SetFromTrueDofs(e0n.Imag());

  // Configure the linear forms for computing S-parameters (projection of the field onto the
  // port mode). Normalize the mode for a chosen polarization direction and unit power,
  // |E x H⋆| ⋅ n, integrated over the port surface (+n is the direction of propagation).
  {
    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    BdrSubmeshHVectorCoefficient<ValueType::REAL> port_nxH0r_func(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0, omega0);
    BdrSubmeshHVectorCoefficient<ValueType::IMAG> port_nxH0i_func(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0, omega0);
    port_sr = std::make_unique<mfem::LinearForm>(&port_nd_fespace->Get());
    port_si = std::make_unique<mfem::LinearForm>(&port_nd_fespace->Get());
    port_sr->AddDomainIntegrator(new VectorFEDomainLFIntegrator(port_nxH0r_func));
    port_si->AddDomainIntegrator(new VectorFEDomainLFIntegrator(port_nxH0i_func));
    port_sr->UseFastAssembly(false);
    port_si->UseFastAssembly(false);
    port_sr->Assemble();
    port_si->Assemble();
    NormalizeWithSign(*port_S0t, *port_E0t, *port_E0n, *port_sr, *port_si);
  }
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeExcitationCoefficientReal() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::REAL>>>(
      attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0,
      omega0);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeExcitationCoefficientImag() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::IMAG>>>(
      attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0,
      omega0);
}

std::unique_ptr<mfem::VectorCoefficient> WavePortData::GetModeFieldCoefficientReal() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshEVectorCoefficient<ValueType::REAL>>>(
      attr_list, *port_E0t, *port_E0n, port_submesh, submesh_parent_elems);
}

std::unique_ptr<mfem::VectorCoefficient> WavePortData::GetModeFieldCoefficientImag() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshEVectorCoefficient<ValueType::IMAG>>>(
      attr_list, *port_E0t, *port_E0n, port_submesh, submesh_parent_elems);
}

double WavePortData::GetExcitationPower() const
{
  // The computed port modes are normalized such that the power integrated over the port is
  // 1: ∫ (E_inc x H_inc⋆) ⋅ n dS = 1.
  return excitation ? 1.0 : 0.0;
}

std::complex<double> WavePortData::GetSParameter(mfem::ParComplexGridFunction &E) const
{
  // Compute port S-parameter, or the projection of the field onto the port mode:
  // (E x H_inc⋆) ⋅ n = E ⋅ (-n x H_inc⋆), integrated over the port surface.
  mfem::ParComplexGridFunction port_E(&port_nd_fespace->Get());
  port_nd_transfer->Transfer(E.real(), port_E.real());
  port_nd_transfer->Transfer(E.imag(), port_E.imag());
  std::complex<double> dot(-((*port_sr) * port_E.real()) - ((*port_si) * port_E.imag()),
                           -((*port_sr) * port_E.imag()) + ((*port_si) * port_E.real()));
  Mpi::GlobalSum(1, &dot, port_nd_fespace->GetComm());
  return dot;
}

std::complex<double> WavePortData::GetPower(mfem::ParComplexGridFunction &E,
                                            mfem::ParComplexGridFunction &B) const
{
  // Compute port power, (E x H) ⋅ n = E ⋅ (-n x H), integrated over the port surface
  // using the computed E and H = μ⁻¹ B fields. The linear form is reconstructed from
  // scratch each time due to changing H. The BdrCurrentVectorCoefficient computes -n x H,
  // where n is an outward normal.
  auto &nd_fespace = *E.ParFESpace();
  const auto &mesh = *nd_fespace.GetParMesh();
  BdrCurrentVectorCoefficient nxHr_func(B.real(), mat_op);
  BdrCurrentVectorCoefficient nxHi_func(B.imag(), mat_op);
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  mfem::LinearForm pr(&nd_fespace), pi(&nd_fespace);
  pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHr_func), attr_marker);
  pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHi_func), attr_marker);
  pr.UseFastAssembly(false);
  pi.UseFastAssembly(false);
  pr.Assemble();
  pi.Assemble();
  std::complex<double> dot(-(pr * E.real()) - (pi * E.imag()),
                           -(pr * E.imag()) + (pi * E.real()));
  Mpi::GlobalSum(1, &dot, nd_fespace.GetComm());
  return dot;
}

WavePortOperator::WavePortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::ParFiniteElementSpace &h1_fespace)
  : suppress_output(false),
    fc(iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, 1.0)),
    kc(1.0 / iodata.DimensionalizeValue(IoData::ValueType::LENGTH, 1.0))
{
  // Set up wave port boundary conditions.
  MFEM_VERIFY(nd_fespace.GetParMesh() == h1_fespace.GetParMesh(),
              "Mesh mismatch in WavePortOperator FE spaces!");
  SetUpBoundaryProperties(iodata, mat_op, nd_fespace, h1_fespace);
  PrintBoundaryInfo(iodata, *nd_fespace.GetParMesh());
}

void WavePortOperator::SetUpBoundaryProperties(const IoData &iodata,
                                               const MaterialOperator &mat_op,
                                               mfem::ParFiniteElementSpace &nd_fespace,
                                               mfem::ParFiniteElementSpace &h1_fespace)
{
  // Check that wave port boundary attributes have been specified correctly.
  const auto &mesh = *nd_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  if (!iodata.boundaries.waveport.empty())
  {
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), port_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    port_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : iodata.boundaries.waveport)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                    "Port boundary attribute tags must be non-negative and correspond to "
                    "boundaries in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1],
                    "Unknown port boundary attribute " << attr << "!");
        MFEM_VERIFY(!port_marker[attr - 1],
                    "Boundary attribute is assigned to more than one wave port!");
        port_marker[attr - 1] = 1;
      }
    }
  }

  // List of all boundaries which will be marked as essential for the purposes of computing
  // wave port modes. This includes all PEC surfaces, but may also include others like when
  // a kinetic inductance or other BC is applied for the 3D simulation but should be
  // considered as PEC for the 2D problem.
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(iodata.boundaries.pec.attributes.size() +
                                   iodata.boundaries.auxpec.attributes.size()));
  for (auto attr : iodata.boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  for (auto attr : iodata.boundaries.auxpec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  // If user accidentally specifies a surface as both "PEC" and "WavePortPEC", this is fine
  // so allow for duplicates in the attribute list.
  dbc_bcs.Sort();
  dbc_bcs.Unique();

  // Set up wave port data structures.
  for (const auto &[idx, data] : iodata.boundaries.waveport)
  {
    ports.try_emplace(idx, data, iodata.solver, mat_op, nd_fespace, h1_fespace, dbc_bcs);
  }
  MFEM_VERIFY(
      ports.empty() || iodata.problem.type == config::ProblemData::Type::DRIVEN,
      "Wave port boundaries are only available for frequency domain driven simulations!");
}

void WavePortOperator::PrintBoundaryInfo(const IoData &iodata, const mfem::ParMesh &mesh)
{
  // Print out BC info for all port attributes.
  if (ports.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin impedance BC for wave ports at attributes:\n");
  for (const auto &[idx, data] : ports)
  {
    for (auto attr : data.GetAttrList())
    {
      const mfem::Vector &normal = data.port_normal;
      Mpi::Print(" {:d}: Index = {:d}, mode = {:d}, d = {:.3e} m", attr, idx, data.mode_idx,
                 iodata.DimensionalizeValue(IoData::ValueType::LENGTH, data.d_offset));
      if (mesh.SpaceDimension() == 3)
      {
        Mpi::Print(", n = ({:+.1f}, {:+.1f}, {:+.1f})", normal(0), normal(1), normal(2));
      }
      else
      {
        Mpi::Print(", n = ({:+.1f}, {:+.1f})", normal(0), normal(1));
      }
      Mpi::Print("\n");
    }
  }

  // Print some information for excited wave ports.
  bool first = true;
  for (const auto &[idx, data] : ports)
  {
    if (!data.excitation)
    {
      continue;
    }
    if (first)
    {
      Mpi::Print("\nConfiguring wave port excitation source term at attributes:\n");
      first = false;
    }
    for (auto attr : data.GetAttrList())
    {
      Mpi::Print(" {:d}: Index = {:d}\n", attr, idx);
    }
  }
}

const WavePortData &WavePortOperator::GetPort(int idx) const
{
  auto it = ports.find(idx);
  MFEM_VERIFY(it != ports.end(), "Unknown wave port index requested!");
  return it->second;
}

mfem::Array<int> WavePortOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    attr_list.Append(data.GetAttrList());
  }
  return attr_list;
}

void WavePortOperator::Initialize(double omega)
{
  bool init = false, first = true;
  for (const auto &[idx, data] : ports)
  {
    init = init || (data.omega0 != omega);
    first = first && (data.omega0 == 0.0);
  }
  if (!init)
  {
    return;
  }
  BlockTimer bt(Timer::WAVEPORT);
  if (!suppress_output)
  {
    Mpi::Print(
        "\nCalculating boundary modes at wave ports for ω/2π = {:.3e} GHz ({:.3e})\n",
        omega * fc, omega);
  }
  for (auto &[idx, data] : ports)
  {
    data.Initialize(omega);
    if (!suppress_output)
    {
      if (first)
      {
        Mpi::Print(" Number of global unknowns for port {:d}:\n"
                   "  H1: {:d}, ND: {:d}\n",
                   idx, data.GlobalTrueH1Size(), data.GlobalTrueNDSize());
      }
      Mpi::Print(" Port {:d}, mode {:d}: kₙ = {:.3e}{:+.3e}i m⁻¹\n", idx, data.mode_idx,
                 data.kn0.real() * kc, data.kn0.imag() * kc);
    }
  }
}

void WavePortOperator::AddExtraSystemBdrCoefficients(double omega,
                                                     MaterialPropertyCoefficient &fbr,
                                                     MaterialPropertyCoefficient &fbi)
{
  // Add wave port boundaries to the bilinear form. This looks a lot like the lumped port
  // boundary, except the iω / Z_s coefficient goes to ikₙ / μ where kₙ is specific to the
  // port mode at the given operating frequency (note only the real part of the propagation
  // constant contributes).
  Initialize(omega);
  for (const auto &[idx, data] : ports)
  {
    const MaterialOperator &mat_op = data.mat_op;
    MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    muinv_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(data.GetAttrList()));
    // fbr.AddCoefficient(muinv_func.GetAttributeToMaterial(),
    //                    muinv_func.GetMaterialProperties(),
    //                    -data.kn0.imag());
    fbi.AddCoefficient(muinv_func.GetAttributeToMaterial(),
                       muinv_func.GetMaterialProperties(), data.kn0.real());
  }
}

void WavePortOperator::AddExcitationBdrCoefficients(double omega, SumVectorCoefficient &fbr,
                                                    SumVectorCoefficient &fbi)
{
  // Re{-U_inc} = Re{+2 (-iω) n x H_inc}, which is a function of E_inc as computed by the
  // modal solution (stored as a grid function and coefficient during initialization).
  // Likewise for the imaginary part.
  Initialize(omega);
  for (const auto &[idx, data] : ports)
  {
    if (!data.excitation)
    {
      continue;
    }
    fbr.AddCoefficient(data.GetModeExcitationCoefficientImag(), 2.0 * omega);
    fbi.AddCoefficient(data.GetModeExcitationCoefficientReal(), -2.0 * omega);
  }
}

}  // namespace palace
