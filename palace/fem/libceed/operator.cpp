// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <numeric>
#include <ceed/backend.h>
#include <mfem/general/forall.hpp>
#include "fem/fespace.hpp"
#include "utils/omp.hpp"

namespace palace::ceed
{

Operator::~Operator()
{
  for (std::size_t i = 0; i < ops.size(); i++)
  {
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(ops[i], &ceed));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&ops[i]));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&ops_t[i]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&u[i]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&v[i]));
  }
}

void Operator::AddOper(CeedOperator op, CeedOperator op_t)
{
  Ceed ceed;
  CeedSize l_in, l_out;
  CeedVector loc_u, loc_v;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op, &ceed));
  PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(op, &l_in, &l_out));
  MFEM_VERIFY((l_in == 0 && l_out == 0) || (mfem::internal::to_int(l_in) == width &&
                                            mfem::internal::to_int(l_out) == height),
              "Dimensions mismatch for CeedOperator!");
  if (op_t)
  {
    CeedSize l_in_t, l_out_t;
    PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(op_t, &l_in_t, &l_out_t));
    MFEM_VERIFY((l_in_t == 0 && l_out_t == 0) || (l_in_t == l_out && l_out_t == l_in),
                "Dimensions mismatch for transpose CeedOperator!");
  }
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, l_in, &loc_u));
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, l_out, &loc_v));

  PalacePragmaOmp(critical(AddOper))
  {
    ops.push_back(op);
    ops_t.push_back(op_t);
    u.push_back(loc_u);
    v.push_back(loc_v);
  }
}

void Operator::AssembleDiagonal(Vector &diag) const
{
  Ceed ceed;
  CeedMemType mem;
  CeedScalar *data;
  MFEM_VERIFY(diag.Size() == height, "Invalid size for diagonal vector!");
  diag = 0.0;
  PalaceCeedCallBackend(CeedOperatorGetCeed(ops[0], &ceed));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    data = diag.ReadWrite();
  }
  else
  {
    data = diag.HostReadWrite();
    mem = CEED_MEM_HOST;
  }

  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < ops.size(); i++)
  {
    Ceed ceed_i;
    PalaceCeedCallBackend(CeedOperatorGetCeed(ops[i], &ceed_i));
    PalaceCeedCall(ceed_i, CeedVectorSetArray(v[i], mem, CEED_USE_POINTER, data));
    PalaceCeedCall(ceed_i, CeedOperatorLinearAssembleAddDiagonal(ops[i], v[i],
                                                                 CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed_i, CeedVectorTakeArray(v[i], mem, nullptr));
  }
}

namespace
{

inline void CeedAddMult(const std::vector<CeedOperator> &ops,
                        const std::vector<CeedVector> &u, const std::vector<CeedVector> &v,
                        const Vector &x, Vector &y)
{
  Ceed ceed;
  CeedMemType mem;
  const CeedScalar *x_data;
  CeedScalar *y_data;
  PalaceCeedCallBackend(CeedOperatorGetCeed(ops[0], &ceed));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    x_data = x.Read();
    y_data = y.ReadWrite();
  }
  else
  {
    x_data = x.HostRead();
    y_data = y.HostReadWrite();
    mem = CEED_MEM_HOST;
  }

  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < ops.size(); i++)
  {
    if (ops[i])  // No-op for an empty operator
    {
      Ceed ceed_i;
      PalaceCeedCallBackend(CeedOperatorGetCeed(ops[i], &ceed_i));
      PalaceCeedCall(ceed_i, CeedVectorSetArray(u[i], mem, CEED_USE_POINTER,
                                                const_cast<CeedScalar *>(x_data)));
      PalaceCeedCall(ceed_i, CeedVectorSetArray(v[i], mem, CEED_USE_POINTER, y_data));
      PalaceCeedCall(ceed_i,
                     CeedOperatorApplyAdd(ops[i], u[i], v[i], CEED_REQUEST_IMMEDIATE));
      PalaceCeedCall(ceed_i, CeedVectorTakeArray(u[i], mem, nullptr));
      PalaceCeedCall(ceed_i, CeedVectorTakeArray(v[i], mem, nullptr));
    }
  }
}

}  // namespace

void Operator::Mult(const Vector &x, Vector &y) const
{
  y = 0.0;
  CeedAddMult(ops, u, v, x, y);
  if (dof_multiplicity.Size() > 0)
  {
    y *= dof_multiplicity;
  }
}

void Operator::AddMult(const Vector &x, Vector &y, const double a) const
{
  MFEM_VERIFY(a == 1.0, "ceed::Operator::AddMult only supports coefficient = 1.0!");
  if (dof_multiplicity.Size() > 0)
  {
    temp.SetSize(height);
    temp = 0.0;
    CeedAddMult(ops, u, v, x, temp);
    if (dof_multiplicity.Size() > 0)
    {
      temp *= dof_multiplicity;
    }
    y += temp;
  }
  else
  {
    CeedAddMult(ops, u, v, x, y);
  }
}

void Operator::MultTranspose(const Vector &x, Vector &y) const
{
  y = 0.0;
  AddMultTranspose(x, y);
}

void Operator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  MFEM_VERIFY(a == 1.0,
              "ceed::Operator::AddMultTranspose only supports coefficient = 1.0!");
  if (dof_multiplicity.Size() > 0)
  {
    temp = x;
    temp *= dof_multiplicity;
    CeedAddMult(ops_t, v, u, temp, y);
  }
  else
  {
    CeedAddMult(ops_t, v, u, x, y);
  }
}

namespace
{

int CeedInternalCallocArray(size_t n, size_t unit, void *p)
{
  *(void **)p = calloc(n, unit);
  MFEM_ASSERT(!n || !unit || *(void **)p,
              "calloc failed to allocate " << n << " members of size " << unit << "!");
  return 0;
}

int CeedInternalFree(void *p)
{
  free(*(void **)p);
  *(void **)p = nullptr;
  return 0;
}

#define CeedInternalCalloc(n, p) CeedInternalCallocArray((n), sizeof(**(p)), p)

void CeedOperatorAssembleCOORemoveZeros(Ceed ceed, CeedSize *nnz, CeedInt **rows,
                                        CeedInt **cols, CeedVector *vals, CeedMemType *mem)
{
  // Filter out zero entries. For now, eliminating zeros happens all on the host.
  // XX TODO: Use Thrust for this (thrust::copy_if and thrust::zip_iterator)
  CeedInt *new_rows, *new_cols;
  PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, &new_rows));
  PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, &new_cols));

  CeedVector new_vals;
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, *nnz, &new_vals));

  CeedSize q = 0;
  const CeedScalar *vals_array;
  CeedScalar *new_vals_array;
  PalaceCeedCall(ceed, CeedVectorGetArrayRead(*vals, CEED_MEM_HOST, &vals_array));
  PalaceCeedCall(ceed, CeedVectorGetArrayWrite(new_vals, CEED_MEM_HOST, &new_vals_array));
  for (CeedSize k = 0; k < *nnz; k++)
  {
    if (vals_array[k] != 0.0)
    {
      new_rows[q] = (*rows)[k];
      new_cols[q] = (*cols)[k];
      new_vals_array[q] = vals_array[k];
      q++;
    }
  }
  PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(*vals, &vals_array));
  PalaceCeedCall(ceed, CeedVectorRestoreArray(new_vals, &new_vals_array));

  PalaceCeedCall(ceed, CeedInternalFree(rows));
  PalaceCeedCall(ceed, CeedInternalFree(cols));
  PalaceCeedCall(ceed, CeedVectorDestroy(vals));

  *rows = new_rows;
  *cols = new_cols;
  *vals = new_vals;
  *nnz = q;
}

void CeedOperatorAssembleCOO(const Operator &op, bool skip_zeros, CeedSize *nnz,
                             CeedInt **rows, CeedInt **cols, CeedVector *vals,
                             CeedMemType *mem)
{
  Ceed ceed;
  CeedScalar *vals_array;
  std::vector<CeedSize> loc_nnz(op.Size()), loc_offsets(op.Size() + 1);
  std::vector<CeedInt *> loc_rows(op.Size()), loc_cols(op.Size());
  std::vector<CeedVector> loc_vals(op.Size());

  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, mem));
  if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) || *mem != CEED_MEM_DEVICE)
  {
    *mem = CEED_MEM_HOST;
  }

  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < op.Size(); i++)
  {
    Ceed ceed_i;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[i], &ceed_i));

    // Assemble sparsity pattern (rows, cols are always host pointers).
    PalaceCeedCall(ceed_i, CeedOperatorLinearAssembleSymbolic(op[i], &loc_nnz[i],
                                                              &loc_rows[i], &loc_cols[i]));

    // Assemble values.
    PalaceCeedCall(ceed_i, CeedVectorCreate(ceed_i, loc_nnz[i], &loc_vals[i]));
    PalaceCeedCall(ceed_i, CeedOperatorLinearAssemble(op[i], loc_vals[i]));
  }

  loc_offsets[0] = 0;
  std::inclusive_scan(loc_nnz.begin(), loc_nnz.end(), loc_offsets.begin() + 1);
  *nnz = loc_offsets.back();
  if (op.Size() == 1)
  {
    // Assemble values.
    *rows = loc_rows[0];
    *cols = loc_cols[0];
    *vals = loc_vals[0];
  }
  else
  {
    // Global assembly.
    PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, rows));
    PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, cols));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, *nnz, vals));
    PalaceCeedCall(ceed, CeedVectorGetArrayWrite(*vals, *mem, &vals_array));

    PalacePragmaOmp(parallel for schedule(static))
    for (std::size_t i = 0; i < op.Size(); i++)
    {
      const auto start = loc_offsets[i];
      const auto end = loc_offsets[i + 1];
      for (auto k = start; k < end; k++)
      {
        (*rows)[k] = loc_rows[i][k - start];
        (*cols)[k] = loc_cols[i][k - start];
      }

      // The CeedVector is on only on device when MFEM is also using the device.
      Ceed ceed_i;
      const CeedScalar *loc_vals_array;
      PalaceCeedCallBackend(CeedVectorGetCeed(loc_vals[i], &ceed_i));
      PalaceCeedCall(ceed_i, CeedVectorGetArrayRead(loc_vals[i], *mem, &loc_vals_array));
      if (*mem != CEED_MEM_HOST)
      {
        mfem::forall(end - start, [=] MFEM_HOST_DEVICE(int k)
                     { vals_array[k + start] = loc_vals_array[k]; });
      }
      else
      {
        for (auto k = start; k < end; k++)
        {
          vals_array[k] = loc_vals_array[k - start];
        }
      }
      PalaceCeedCall(ceed_i, CeedVectorRestoreArrayRead(loc_vals[i], &loc_vals_array));
      PalaceCeedCall(ceed_i, CeedInternalFree(&loc_rows[i]));
      PalaceCeedCall(ceed_i, CeedInternalFree(&loc_cols[i]));
      PalaceCeedCall(ceed_i, CeedVectorDestroy(&loc_vals[i]));
    }

    PalaceCeedCall(ceed, CeedVectorRestoreArray(*vals, &vals_array));
  }

  // std::cout << "  Operator full assembly (COO) has " << *nnz << " NNZ";
  if (skip_zeros && *nnz > 0)
  {
    CeedOperatorAssembleCOORemoveZeros(ceed, nnz, rows, cols, vals, mem);
    // std::cout << " (new NNZ after removal: " << *nnz << ")";
  }
  // std::cout << "\n";
}

}  // namespace

std::unique_ptr<mfem::SparseMatrix> CeedOperatorFullAssemble(const Operator &op,
                                                             bool skip_zeros, bool set)
{
  // First, get matrix on master thread in COO format, withs rows/cols always on host and
  // vals potentially on the device. Process skipping zeros if desired.
  Ceed ceed;
  CeedSize nnz;
  CeedInt *rows, *cols;
  CeedVector vals;
  CeedMemType mem;
  CeedOperatorAssembleCOO(op, skip_zeros, &nnz, &rows, &cols, &vals, &mem);
  PalaceCeedCallBackend(CeedVectorGetCeed(vals, &ceed));

  // Preallocate CSR memory on host (like PETSc's MatSetValuesCOO).
  auto mat = std::make_unique<mfem::SparseMatrix>();
  mat->OverrideSize(op.Height(), op.Width());
  mat->GetMemoryI().New(op.Height() + 1);
  auto *I = mat->GetI();
  mfem::Array<int> J(nnz), perm(nnz), Jmap(nnz + 1);

  for (int i = 0; i < op.Height() + 1; i++)
  {
    I[i] = 0;
  }
  for (int k = 0; k < nnz; k++)
  {
    perm[k] = k;
  }
  std::sort(perm.begin(), perm.end(),
            [&](const int &i, const int &j) { return (rows[i] < rows[j]); });

  int q = -1;  // True nnz index
  for (int k = 0; k < nnz;)
  {
    // Sort column entries in the row.
    const int row = rows[perm[k]];
    const int start = k;
    while (k < nnz && rows[perm[k]] == row)
    {
      k++;
    }
    std::sort(perm.begin() + start, perm.begin() + k,
              [&](const int &i, const int &j) { return (cols[i] < cols[j]); });

    q++;
    I[row + 1] = 1;
    J[q] = cols[perm[start]];
    Jmap[q + 1] = 1;
    for (int p = start + 1; p < k; p++)
    {
      if (cols[perm[p]] != cols[perm[p - 1]])
      {
        // New nonzero.
        q++;
        I[row + 1]++;
        J[q] = cols[perm[p]];
        Jmap[q + 1] = 1;
      }
      else
      {
        Jmap[q + 1]++;
      }
    }
  }
  PalaceCeedCall(ceed, CeedInternalFree(&rows));
  PalaceCeedCall(ceed, CeedInternalFree(&cols));
  const int nnz_new = q + 1;

  // Finalize I, Jmap.
  I[0] = 0;
  for (int i = 0; i < op.Height(); i++)
  {
    I[i + 1] += I[i];
  }
  Jmap[0] = 0;
  for (int k = 0; k < nnz_new; k++)
  {
    Jmap[k + 1] += Jmap[k];
  }

  mat->GetMemoryJ().New(nnz_new, mat->GetMemoryJ().GetMemoryType());
  mat->GetMemoryData().New(nnz_new, mat->GetMemoryJ().GetMemoryType());
  {
    const auto *d_J_old = J.Read();
    auto *d_J = mfem::Write(mat->GetMemoryJ(), nnz_new);
    mfem::forall(nnz_new, [=] MFEM_HOST_DEVICE(int k) { d_J[k] = d_J_old[k]; });
  }

  // Fill the values (on device).
  const CeedScalar *vals_array;
  PalaceCeedCall(ceed, CeedVectorGetArrayRead(vals, mem, &vals_array));
  {
    const auto *d_perm = perm.Read();
    const auto *d_Jmap = Jmap.Read();
    auto *d_A = mfem::Write(mat->GetMemoryData(), nnz_new);
    if (set)
    {
      mfem::forall(nnz_new,
                   [=] MFEM_HOST_DEVICE(int k) { d_A[k] = vals_array[d_perm[d_Jmap[k]]]; });
    }
    else
    {
      mfem::forall(nnz_new,
                   [=] MFEM_HOST_DEVICE(int k)
                   {
                     double sum = 0.0;
                     for (int p = d_Jmap[k]; p < d_Jmap[k + 1]; p++)
                     {
                       sum += vals_array[d_perm[p]];
                     }
                     d_A[k] = sum;
                   });
    }
  }
  PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(vals, &vals_array));
  PalaceCeedCall(ceed, CeedVectorDestroy(&vals));

  return mat;
}

std::unique_ptr<Operator> CeedOperatorCoarsen(const Operator &op_fine,
                                              const FiniteElementSpace &fespace_coarse)
{
  auto SingleOperatorCoarsen =
      [&fespace_coarse](Ceed ceed, CeedOperator op_fine, CeedOperator *op_coarse)
  {
    CeedBasis basis_fine;
    CeedElemTopology geom;
    PalaceCeedCall(ceed, CeedOperatorGetActiveBasis(op_fine, &basis_fine));
    PalaceCeedCall(ceed, CeedBasisGetTopology(basis_fine, &geom));

    const auto &geom_data =
        fespace_coarse.GetMesh().GetCeedGeomFactorData(ceed).at(GetMfemTopology(geom));
    CeedElemRestriction restr_coarse = fespace_coarse.GetCeedElemRestriction(
        ceed, GetMfemTopology(geom), geom_data.indices);
    CeedBasis basis_coarse = fespace_coarse.GetCeedBasis(ceed, GetMfemTopology(geom));

    PalaceCeedCall(ceed, CeedOperatorMultigridLevelCreate(op_fine, nullptr, restr_coarse,
                                                          basis_coarse, op_coarse, nullptr,
                                                          nullptr));
  };

  // Initialize the coarse operator.
  auto op_coarse = std::make_unique<SymmetricOperator>(fespace_coarse.GetVSize(),
                                                       fespace_coarse.GetVSize());

  // Assemble the coarse operator by coarsening each sub-operator (over threads, geometry
  // types, integrators) of the original fine operator.
  MFEM_VERIFY(internal::GetCeedObjects().size() == op_fine.Size(),
              "Unexpected size mismatch in multithreaded Ceed contexts!");
  const std::size_t nt = op_fine.Size();
  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < nt; i++)
  {
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op_fine[i], &ceed));
    {
      Ceed ceed_parent;
      PalaceCeedCall(ceed, CeedGetParent(ceed, &ceed_parent));
      if (ceed_parent)
      {
        ceed = ceed_parent;
      }
    }

    // Initialize the composite operator on each thread.
    CeedOperator loc_op;
    PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op));

    bool composite;
    PalaceCeedCall(ceed, CeedOperatorIsComposite(op_fine[i], &composite));
    if (composite)
    {
      CeedInt nloc_ops_fine;
      CeedOperator *loc_ops_fine;
      PalaceCeedCall(ceed, CeedCompositeOperatorGetNumSub(op_fine[i], &nloc_ops_fine));
      PalaceCeedCall(ceed, CeedCompositeOperatorGetSubList(op_fine[i], &loc_ops_fine));
      for (CeedInt k = 0; k < nloc_ops_fine; k++)
      {
        CeedOperator sub_op;
        SingleOperatorCoarsen(ceed, loc_ops_fine[k], &sub_op);
        PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
        PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
      }
    }
    else
    {
      CeedOperator sub_op;
      SingleOperatorCoarsen(ceed, op_fine[i], &sub_op);
      PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
      PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
    }
    PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op));
    op_coarse->AddOper(loc_op);  // Thread-safe
  }

  return op_coarse;
}

}  // namespace palace::ceed
