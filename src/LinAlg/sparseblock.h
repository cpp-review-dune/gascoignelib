/*----------------------------   sparseblock.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __sparseblock_H
#define __sparseblock_H
/*----------------------------   sparseblock.h     ---------------------------*/

#include "../Common/entrymatrix.h"
#include "../Interface/gascoigne.h"

#include "sparse_20.h"

/*****************************************************/

namespace Gascoigne {
class SparseBlock : public numfixarray<SPARSE_NENTRIES, MatrixEntryType>
{
  typedef nvector<double>::iterator viterator;
  typedef nvector<double>::const_iterator const_viterator;

public:
  ShortIndexType ncomp() const { return SPARSE_NROWS; }

  virtual void operator*=(const SparseBlock&) = 0;
  virtual void operator*=(double s) = 0;

  virtual void transpose() = 0;
  virtual void transpose(SparseBlock& A) = 0;
  virtual void copy_transpose(const SparseBlock& A) = 0;

  virtual void zero_row(IndexType) = 0;
  virtual void uno_diag(IndexType) = 0;
  virtual MatrixEntryType& diag(IndexType i) = 0;
  virtual void getrow(std::vector<double>& v, IndexType i) = 0;
  virtual void getcolumn(std::vector<double>& v, IndexType i) = 0;
  virtual void setrow(std::vector<double>& v, IndexType i) = 0;
  virtual void setcolumn(std::vector<double>& v, IndexType i) = 0;

  double operator()(IndexType r, IndexType c) const
  {
    IndexType p = SPARSE_START[r];
    for (; p < SPARSE_START[r + 1]; ++p)
      if (SPARSE_COL[p] == c)
        break;
    if (p == SPARSE_START[r + 1])
      assert(0);
    return 0;
    return numfixarray<SPARSE_NENTRIES, MatrixEntryType>::operator[](p);
  }

  void DirichletRow(const std::vector<IndexType>& cv)
  {
    for (auto row : cv) {
      for (IndexType p = SPARSE_START[row]; p < SPARSE_START[row + 1]; ++p)
        numfixarray<SPARSE_NENTRIES, MatrixEntryType>::operator[](p) = 0;
    }
  }

  void DirichletCol(const std::vector<IndexType>& cv)
  {
    for (auto col : cv) {
      for (IndexType r = 0; r < SPARSE_NROWS; ++r)
        for (IndexType p = SPARSE_START[r]; p < SPARSE_START[r + 1]; ++p)
          if (SPARSE_COL[p] == col)
            numfixarray<SPARSE_NENTRIES, MatrixEntryType>::operator[](p) = 0;
    }
  }

  void DirichletDiag(const std::vector<IndexType>& cv)
  {
    for (auto row : cv) {
      IndexType p;
      for (p = SPARSE_START[row]; p < SPARSE_START[row + 1]; ++p)
        if (SPARSE_COL[p] == row)
          break;

      assert(p < SPARSE_START[row + 1]);
      numfixarray<SPARSE_NENTRIES, MatrixEntryType>::operator[](p) = 1.;
    }
  }

  virtual void entry(const nmatrix<double>&) = 0;

  void entry(IndexType i, IndexType j, const EntryMatrix& E, double s = 1.)
  {
    for (IndexType row = 0; row < SPARSE_NROWS; ++row)
      for (IndexType p = SPARSE_START[row]; p < SPARSE_START[row + 1]; ++p) {
        IndexType col = SPARSE_COL[p];
#pragma omp atomic update
        numfixarray<SPARSE_NENTRIES, MatrixEntryType>::operator[](p) +=
          s * E(i, j, row, col);
      }
  }

  void dual_entry(IndexType i, IndexType j, const EntryMatrix&, double s = 1.)
  {
    assert(0);
  }
  virtual void inverse() = 0;
  virtual void vmult(viterator) const = 0;
  virtual void mult(SparseBlock&, const SparseBlock&) const = 0;

  virtual void submult(const SparseBlock& B, const SparseBlock& C) = 0;

  virtual void add(double s, const SparseBlock& A) = 0;
  virtual void adddiag(const nvector<double>& s, double l) = 0;

  virtual void add(double s, const TimePattern& TP) = 0;

  void cadd(double s, viterator p, const_viterator q0) const
  {
    auto it = this->begin();
    //	numfixarray<STRUCT::NENTRIES,MatrixEntryType>::begin();
    for (IndexType row = 0; row < SPARSE_NROWS; ++row) {
      double sum = 0.0;

      for (IndexType pos = SPARSE_START[row]; pos < SPARSE_START[row + 1];
           ++pos, ++it) {
        IndexType col = SPARSE_COL[pos];
        sum += (*it) * *(q0 + col);
      }

      *p += s * sum;
      p++;
    }
  }

  virtual void caddtrans(double s, viterator p, const_viterator q0) const = 0;
  virtual void subtract(viterator p0, const_viterator q0) const = 0;
  virtual std::ostream& print(std::ostream& s) const = 0;

  // Zugriff auf Inhalt ueber ganzen Vektor, damits auch ohne
  // Struktur geht.
  virtual void vector_get(nvector<MatrixEntryType>& v) const = 0;
  virtual void vector_set(nvector<MatrixEntryType>& v) = 0;
  virtual void vector_add(double d, nvector<MatrixEntryType>& v) = 0;

  friend std::ostream& operator<<(std::ostream& s, const SparseBlock& A)
  {
    return s;
  }
};

} // namespace Gascoigne

/*----------------------------   sparseblock.h     ---------------------------*/
/* end of #ifndef __sparseblock_H */
#endif
/*----------------------------   sparseblock.h     ---------------------------*/
