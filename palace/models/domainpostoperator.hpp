// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
#define PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP

#include <map>
#include <memory>
#include <utility>
#include <mfem.hpp>
#include "linalg/operator.hpp"

namespace palace
{

class IoData;
class MaterialOperator;

//
// A class handling domain postprocessing.
//
class DomainPostOperator
{
private:
  // Bilinear forms for computing field energy integrals over domains.
  std::unique_ptr<Operator> M_ND, M_RT;
  std::map<int, std::pair<std::unique_ptr<Operator>, std::unique_ptr<Operator>>> M_i;

  // Temporary vectors for inner product calculations.
  mutable Vector D, H;

public:
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const mfem::ParFiniteElementSpace *nd_fespace,
                     const mfem::ParFiniteElementSpace *rt_fespace);

  // Access data structures for the postprocessing the domain with the given type.
  const auto &GetBulkDomains() const { return M_i; }

  // Get volume integrals computing the electric or magnetic field energy in the domain.
  double GetElectricFieldEnergy(const mfem::ParComplexGridFunction &E) const;
  double GetElectricFieldEnergy(const mfem::ParGridFunction &E) const;
  double GetMagneticFieldEnergy(const mfem::ParComplexGridFunction &B) const;
  double GetMagneticFieldEnergy(const mfem::ParGridFunction &B) const;

  // Get volume integrals for bulk electric or magnetic field energy in a portion of the
  // domain.
  double GetDomainElectricFieldEnergy(int idx, const mfem::ParComplexGridFunction &E) const;
  double GetDomainElectricFieldEnergy(int idx, const mfem::ParGridFunction &E) const;
  double GetDomainMagneticFieldEnergy(int idx, const mfem::ParComplexGridFunction &E) const;
  double GetDomainMagneticFieldEnergy(int idx, const mfem::ParGridFunction &E) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
