/**
*
* Copyright (C) 2004, 2009, 2011 by the Gascoigne 3D authors
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


#ifndef  __PointMatrix_h
#define  __PointMatrix_h

#include  "matrixinterface.h"
#include  "simplematrix.h"
#include  "sparsestructureadaptor.h"
#include  "mginterpolatormatrix.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments PointMatrix

///
///
/////////////////////////////////////////////

class PointMatrix : public SimpleMatrix, virtual public MatrixInterface
{
private:
  template<bool atom>
  void entry_universal(nvector<int>::const_iterator start1,
                       nvector<int>::const_iterator stop1,
                       nvector<int>::const_iterator start2,
                       nvector<int>::const_iterator stop2, const EntryMatrix& M,
                       double s = 1.);
  template<bool atom>
  void entry_universal(niiterator start, niiterator stop, const EntryMatrix& M,
                       double s = 1.);

protected:

  int _ncomp;
  SparseStructureAdaptor* SSAP;

public:

//
///  Constructor
//
    PointMatrix(int ncomp, std::string type);
    virtual ~PointMatrix();

    std::string GetName() const {return "PointMatrix";}

    void zero() {
      SimpleMatrix::zero();
    }
    void vmult(GlobalVector& y, const GlobalVector& x, double d=1.) const;
    void vmult_transpose(GlobalVector& y, const GlobalVector& x, double d=1.) const;

    const StencilInterface* GetStencil() const { return SimpleMatrix::GetStencil();}
    void ReInit(const SparseStructureInterface* S);

    void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.);
  void entry(niiterator start1, niiterator stop1,
	     niiterator start2, niiterator stop2,
	     const EntryMatrix& M, double s=1.);
  void entry_atomic(nvector<int>::const_iterator start1,
                    nvector<int>::const_iterator stop1,
                    nvector<int>::const_iterator start2,
                    nvector<int>::const_iterator stop2, const EntryMatrix& M,
                    double s = 1.);
  void entry_atomic(niiterator start, niiterator stop, const EntryMatrix& M,
                    double s = 1.);
  void entry_diag(int i, const nmatrix<double>& M);
  void dirichlet(int i, const std::vector<int>& cv);
  void dirichlet_only_row(int i, const std::vector<int>& cv);
  void periodic(const std::map<int, int>& m_PeriodicPairs,
                const IntVector& iv_Components);

  void transpose()
  {
    SimpleMatrix::transpose();
    }

    void AddMassWithDifferentStencil(const MatrixInterface* M, const TimePattern& TP, double s=1.);
    void AddMassWithDifferentStencilJacobi(const MatrixInterface* M, const TimePattern& TP, double s=1.);

    void RestrictMatrix(const MgInterpolatorMatrix& I, const PointMatrix& Ah);
};
  template<bool atom>
  void PointMatrix::entry_universal(niiterator start,
                          niiterator stop,
                          const EntryMatrix &M,
                          double s)
  {
    int n = stop - start;

    for (int ii = 0; ii < n; ii++)
    {
      int i = *(start + ii);
      for (int c = 0; c < _ncomp; c++)
      {
        int iglob = SSAP->index(i, c);
        for (int jj = 0; jj < n; jj++)
        {
          int j = *(start + jj);
          for (int d = 0; d < _ncomp; d++)
          {
            int jglob = SSAP->index(j, d);
            int pos = ST.Find(iglob, jglob);
            //if constexpr (atom)
            if (atom)
            {
#pragma omp atomic update
              value[pos] += s * M(ii, jj, c, d);
            }
            else
            {
              value[pos] += s * M(ii, jj, c, d);
            }
          }
        }
      }
    }
  }

  template<bool atom>
  void PointMatrix::entry_universal(niiterator start1,
                          niiterator stop1,
                          niiterator start2,
                          niiterator stop2,
                          const EntryMatrix &M,
                          double s)
  {
    int n1 = stop1 - start1;
    int n2 = stop2 - start2;
    assert(n1==n2);


    for (int ii = 0; ii < n1; ii++)
    {
      int i = *(start1 + ii);
      for (int c = 0; c < _ncomp; c++)
      {
        int iglob = SSAP->index(i, c);
        for (int jj = 0; jj < n1; jj++)
        {
          int j = *(start2 + jj);
          for (int d = 0; d < _ncomp; d++)
          {
            int jglob = SSAP->index(j, d);
            int pos = ST.Find(iglob, jglob);
            //if constexpr(atom)
            if(atom)
            {
#pragma omp atomic update
              value[pos] += s * M(ii, jj, c, d);
            }
            else
            {
              value[pos] += s * M(ii, jj, c, d);
            }
          }
        }
      }
    }
  }
}

#endif
