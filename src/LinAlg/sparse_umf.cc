#include "sparse_umf.xx"
#include "fmatrixblock.h"
#include "sparseblock.h"

#ifdef __WITH_UMFPACK__

namespace Gascoigne {

template<>
void
SparseUmf<SparseBlock>::ConstructStructure(const IntVector& perm,
                                           const MatrixInterface& __XXX)
{
  assert(__AS);
  // reinit size of UMFPACK-structure
  const ColumnStencil* ST =
    dynamic_cast<const ColumnStencil*>(__AS->GetStencil());
  assert(static_cast<const void*>(__AS) == static_cast<const void*>(&__XXX));

  assert(ST);
  __ncomp = __AS->mat(0)->ncomp();
  assert(__ncomp == SPARSE_NROWS);

  int n = ST->n();
  int nentries = ST->nentries();

  // reserve size for matrix-entries
  __Ax.resize(SPARSE_NENTRIES * nentries);
  __Ax.zero();
  __Ac.clear();
  __Ap.clear();
  __Ap.push_back(0);

  for (int rA = 0; rA < n; ++rA)
  // insert myblock at r,c
  // row-wise, sorted in row
  {
    // use structure indicated in LDBLOCKSTRUCT
    for (int rb = 0; rb < SPARSE_NROWS; ++rb) {
      for (int pA = ST->start(rA); pA != ST->stop(rA); ++pA) {
        int cA = ST->col(pA);
        for (int p = SPARSE_START[rb]; p < SPARSE_START[rb + 1]; ++p)
          __Ac.push_back(__ncomp * cA + SPARSE_COL[p]);
      }
      __Ap.push_back(__Ac.size());
    }
  }
  long n_umf = n * __ncomp;
  long n_entries_umf = SPARSE_NENTRIES * nentries;

  assert(__Ac.size() == n_entries_umf);
  assert(__Ax.size() == n_entries_umf);
  assert(__Ap.size() == n_umf + 1);
  assert(n * SPARSE_NROWS + 1 == __Ap.size());
  assert(__Ap.size() > 0);
  assert(__Ap[__Ap.size() - 1] == n_entries_umf);

  umfpack_dl_free_symbolic(&Symbolic);

  const long* sb = &(*__Ap.begin());
  const long* cb = &(*__Ac.begin());

  int status =
    umfpack_dl_symbolic(n_umf, n_umf, sb, cb, NULL, &Symbolic, Control, Info);

  if (status != UMFPACK_OK) {
    umfpack_dl_report_info(Control, Info);
    umfpack_dl_report_status(Control, status);
    cerr << "umfpack_symbolic failed\n";
    exit(1);
  }
}

template<>
void
SparseUmf<SparseBlock>::copy_entries(const MatrixInterface& A)
{

  // check if __AS and A are the same object...
  // then, we do not need __AS!!!
  assert(static_cast<const void*>(&__AS) == static_cast<const void*>(&A));

  // Copy Entries
  assert(__ncomp == __AS->mat(0)->ncomp());

  const ColumnStencil* ST =
    dynamic_cast<const ColumnStencil*>(__AS->GetStencil());
  assert(ST);
  int pp = 0;
  for (int rA = 0; rA < ST->n(); ++rA)
    for (int rn = 0; rn < SPARSE_NROWS; ++rn) {
      for (int pA = ST->start(rA); pA != ST->stop(rA); ++pA) {
        const SparseBlock& b = *__AS->mat(pA);
        for (int p = SPARSE_START[rn]; p < SPARSE_START[rn + 1]; ++p, ++pp)
          __Ax[pp] = b[p];
      }
    }
  assert(pp == __Ax.size());
}

template class SparseUmf<FMatrixBlock<1>>;
template class SparseUmf<FMatrixBlock<2>>;
template class SparseUmf<FMatrixBlock<3>>;
template class SparseUmf<FMatrixBlock<4>>;
template class SparseUmf<FMatrixBlock<5>>;
template class SparseUmf<FMatrixBlock<6>>;
template class SparseUmf<FMatrixBlock<7>>;

template class SparseUmf<SparseBlock>;

} // namespace Gascoigne

#endif
