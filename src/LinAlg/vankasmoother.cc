#include "vankasmoother.h"
#include "gascoignehash.h"

namespace Gascoigne
{
//////////////////// Construction
void VankaSmoother::ConstructStructure(const IntVector& perm, const MatrixInterface& A)
{
  assert(_dofhandler);
  const SparseBlockMatrix<FMatrixBlock<1>>* M1 =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<1>>*>(&A);
  const SparseBlockMatrix<FMatrixBlock<2>>* M2 =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<2>>*>(&A);
  const SparseBlockMatrix<FMatrixBlock<3>>* M3 =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<3>>*>(&A);
  const SparseBlockMatrix<FMatrixBlock<4>>* M4 =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<4>>*>(&A);
  const SparseBlockMatrix<FMatrixBlock<5>>* M5 =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<5>>*>(&A);
  const SparseBlockMatrix<FMatrixBlock<6>>* M6 =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<6>>*>(&A);
  const SparseBlockMatrix<FMatrixBlock<7>>* M7 =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<7>>*>(&A);
  if (!(M1 || M2 || M3 || M4 || M5 || M6 || M7))
  {
    std::cerr << "VankaSmoother: wrong matrix" << std::endl;
    abort();
  }

  if (M1)
    _ncomp = 1;
  else if (M2)
    _ncomp = 2;
  else if (M3)
    _ncomp = 3;
  else if (M4)
    _ncomp = 4;
  else if (M5)
    _ncomp = 5;
  else if (M6)
    _ncomp = 6;
  else if (M7)
    _ncomp = 7;
  else
    assert(0);

  ////////// Construct patches
  int patchlevel = 2;

  int npatches = _dofhandler->nelements(patchlevel);
  _sizeofpatch = _dofhandler->nodes_per_element(patchlevel);
  _patchlist.resize(npatches, std::vector<int>(_sizeofpatch));
#pragma omp parallel for schedule(static)
  for (int p = 0; p < npatches; ++p)
  {
    const IntVector iop = _dofhandler->GetElement(patchlevel, p);
    assert(iop.size() == _patchlist[p].size());
    for (int i = 0; i < iop.size(); ++i)
      _patchlist[p][i] = iop[i];
  }

  /////////// Init LU-List
  _lu.resize(npatches, Eigen::PartialPivLU<Eigen::MatrixXd>(_sizeofpatch * _ncomp));

  ElementColoring(patchlevel);
}

template <int NCOMP>
void VankaSmoother::copy_entries_sparseblockmatrix(
  const SparseBlockMatrix<FMatrixBlock<NCOMP>>& A)
{
  const ColumnDiagStencil& S = dynamic_cast<const ColumnDiagStencil&>(*A.GetStencil());

  // Copy entries & assemble ILU
#pragma omp parallel for schedule(static)
  for (int p = 0; p < _patchlist.size(); ++p)
  {
    assert(_patchlist[p].size() == _sizeofpatch);

    Eigen::MatrixXd Matrix_on_Block;
    Matrix_on_Block.resize(_sizeofpatch * NCOMP, _sizeofpatch * NCOMP);
    Matrix_on_Block.setZero();

    // inverse index set? store globally?
    HASHMAP<int, int> INP;
    for (int i = 0; i < _patchlist[p].size(); ++i)
      INP[_patchlist[p][i]] = i;

    for (int r = 0; r < _sizeofpatch; ++r)
    {
      assert(r < _patchlist[p].size());
      int row = _patchlist[p][r];
      assert(row < S.n());

      // Copy Matrix
      for (int pos = S.start(row); pos < S.stop(row); ++pos)
      {
        int col = S.col(pos);
        /* columnsinpatch.insert(col); */

        // find inverse index. Skip, if not n patch
        auto inverseindex = INP.find(col);
        if (inverseindex == INP.end())
          continue;
        int c = inverseindex->second;

        const FMatrixBlock<NCOMP>& B = (*A.mat(pos));
        for (int cr = 0; cr < NCOMP; ++cr)
          for (int cc = 0; cc < NCOMP; ++cc)
            Matrix_on_Block(NCOMP * r + cr, NCOMP * c + cc) = B(cr, cc);
      }
    }
    // Compute LU
    _lu[p].compute(Matrix_on_Block);
  }
}

void VankaSmoother::copy_entries(const MatrixInterface& A)
{
  if (_ncomp == 1)
    copy_entries_sparseblockmatrix(
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<1>>&>(A));
  else if (_ncomp == 2)
    copy_entries_sparseblockmatrix(
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<2>>&>(A));
  else if (_ncomp == 3)
    copy_entries_sparseblockmatrix(
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<3>>&>(A));
  else if (_ncomp == 4)
    copy_entries_sparseblockmatrix(
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<4>>&>(A));
  else if (_ncomp == 5)
    copy_entries_sparseblockmatrix(
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<5>>&>(A));
  else if (_ncomp == 6)
    copy_entries_sparseblockmatrix(
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<6>>&>(A));
  else if (_ncomp == 7)
    copy_entries_sparseblockmatrix(
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<7>>&>(A));
}

void VankaSmoother::solve(GlobalVector& x) const
{
  assert(x.ncomp() == _ncomp);
  ////////////////////// Jacobi smoother

  // copy original vector
  GlobalVector B = x;
  // set solution to zero
  x.zero();
  // vector for averaging
  std::vector<int> count(x.n(), 0);

  // perform Vanka loop
//#pragma omp parallel for
// for (int p=0;p<_patchlist.size();++p)
//#pragma omp parallel for
// for (int p=0;p<_patchlist.size();++p)
#pragma omp parallel
  {
    for (int col = 0; col < NumberofColors(); col++)
    {
      const std::vector<int>& ewcol = elementswithcolor(col);
#pragma omp for schedule(static)
      for (int iii = 0; iii < ewcol.size(); ++iii)
      {
        int p = ewcol[iii];
        // copy local patch-vector
        Eigen::VectorXd H(_sizeofpatch * _ncomp);
        for (int r = 0; r < _sizeofpatch; ++r)
          for (int c = 0; c < _ncomp; ++c)
            H(_ncomp * r + c, 0) = B(_patchlist[p][r], c);

        // perform inversion
        H = _lu[p].permutationP() * H;
        _lu[p].matrixLU().triangularView<Eigen::UnitLower>().solveInPlace(H);
        _lu[p].matrixLU().triangularView<Eigen::Upper>().solveInPlace(H);

        // update
        for (int r = 0; r < _sizeofpatch; ++r)
        {
          //#pragma omp critical
          {
            count[_patchlist[p][r]]++;
            for (int c = 0; c < _ncomp; ++c)
              x(_patchlist[p][r], c) += H(_ncomp * r + c, 0);
          }
        }
      }
    }

    // average
#pragma omp for schedule(static)
    for (int i = 0; i < x.n(); ++i)
    {
      assert(count[i] > 0);
      x.scale_node(i, 1.0 / count[i]);
    }
  }
}

}  // namespace Gascoigne
