/**
 *
 * Copyright (C) 2004, 2005, 2006, 2009, 2011 by the Gascoigne 3D authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/

#ifndef __sparseblockmatrix_h
#define __sparseblockmatrix_h

#include "columndiagstencil.h"
#include "gascoigne.h"
#include "matrixinterface.h"
#include "sparsestructure.h"

/*-------------------------------------------------------------*/

namespace Gascoigne {
template<class B>
class SparseBlockMatrix : public virtual MatrixInterface
{
private:
protected:
  typedef typename std::vector<B>::const_iterator const_iterator;
  typedef typename std::vector<B>::iterator iterator;
  typedef typename std::pair<int, int> IntPair;

  ColumnDiagStencil US;
  std::vector<B> smat;
  int nc;

  void matrix_vector_trans(int p,
                           double* yp,
                           const double* xp,
                           double s = 1.) const;

public:
  int size() const { return smat.size(); }

  SparseBlockMatrix();
  SparseBlockMatrix(const SparseBlockMatrix<B>& A);
  virtual ~SparseBlockMatrix() {}

  void transpose();

  std::string GetName() const { return "SparseBlockMatrix"; }

  /////// Zugriff //////////////////////

  const_iterator mat(int pos) const
  {
    assert(pos < smat.size());
    return smat.begin() + pos;
  }
  iterator mat(int pos)
  {
    assert(pos < smat.size());
    return smat.begin() + pos;
  }

  const StencilInterface* GetStencil() const { return &US; }

  int n() const { return US.n(); };
  int nentries() const { return US.nentries(); };
  int ntotal() const { return smat.size(); };

  int rowsize(int i) const { return US.start(i + 1) - US.start(i); }
  const std::vector<B>& mat() const { return smat; }

  ///// Methods //////////////////////

  void AddMassWithDifferentStencil(const MatrixInterface* M,
                                   const TimePattern& TP,
                                   double s = 1.);
  void AddMassWithDifferentStencilJacobi(const MatrixInterface* M,
                                         const TimePattern& TP,
                                         double s = 1.);

  void copy_entries(const MatrixInterface& S);

  SparseBlockMatrix& operator=(const SparseBlockMatrix<B>& S);

  void ReInit(const SparseStructureInterface*);
  void scale_diag(int i, const std::vector<int>& cv, double s);
  void dirichlet(int i, const std::vector<int>& cv);
  void dirichlet_only_row(int i, const std::vector<int>& cv);
  void dirichlet_only_column(int i, const std::vector<int>& cv);
  void dirichlet_only_row_no_diag(int i, const std::vector<int>& cv);
  void periodic(const std::map<int, int>& m_PeriodicPairs,
                const IntVector& iv_Components);

  void zero();
  void entry(nvector<int>::const_iterator start1,
             nvector<int>::const_iterator stop1,
             nvector<int>::const_iterator start2,
             nvector<int>::const_iterator stop2,
             const EntryMatrix& M,
             double s = 1.);
  void entry(nvector<int>::const_iterator start,
             nvector<int>::const_iterator stop,
             const EntryMatrix& M,
             double s = 1.);
  void entrydual(nvector<int>::const_iterator start,
                 nvector<int>::const_iterator stop,
                 const EntryMatrix& M,
                 double s = 1.);

  void GaussSeidel(GlobalVector& y, const GlobalVector& x) const;
  void Jacobi(GlobalVector& x) const;

  void MatrixResidualSome(const std::vector<int>& indices,
                          GlobalVector& h,
                          const GlobalVector& x,
                          const GlobalVector& y) const;
  void vmult(GlobalVector& y, const GlobalVector& x, double s = 1.) const;
  void vmult(GlobalVector& y,
             const GlobalVector& x,
             const TimePattern& TP,
             double s = 1.) const;
  void entry_diag(int i, const nmatrix<double>& M);

  /*-----------------------------------------------*/

  void FillInterfaceList(const nvector<int>& elements,
                         nvector<int>& start,
                         nvector<MatrixEntryType>& values) const;
  void FurbishInterface(double d,
                        const nvector<int>& elements,
                        const nvector<int>& start,
                        const nvector<MatrixEntryType>& values);

  /*-----------------------------------------------*/

  std::ostream& Write(std::ostream& s) const;

  /// write matrix to file in simple format:
  /// row col value
  /// can be used by matlab
  ///   load mat.dat
  ///   A = spconvert(mat)
  void write_raw(std::string fname) const;

  friend std::ostream& operator<<(std::ostream& s,
                                  const SparseBlockMatrix<B>& A)
  {
    std::cerr << "\"ostream& operator<<(ostream &s, const "
                 "SparseBlockMatrix<B>& A)\" not written!"
              << std::endl;
    abort();
  }
};
} // namespace Gascoigne

#endif
