/**
 *
 * Copyright (C) 2004, 2005, 2010 by the Gascoigne 3D authors
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

#ifndef __FiniteElement_h
#define __FiniteElement_h

#include "baseq12d.h"
#include "baseq13d.h"
#include "baseq22d.h"
#include "baseq23d.h"
#include "feminterface.h"
#include "transformation2d.h"
#include "transformation3d.h"

/*-----------------------------------------------------*/

namespace Gascoigne {

/////////////////////////////////////////////
///
///@brief
///  FE based on Transformation (TRAFO) and Referenzelement (BASE)

///
///
/////////////////////////////////////////////

template<int DIM, int BDIM, class TRAFO, class BASE>
class FiniteElement : public FemInterface
{
protected:
  TRAFO T;
  BASE B;
  mutable std::vector<Vertex<DIM>> grad;
  mutable double det;

  virtual void ComputeGrad() const;

public:
  FiniteElement();

  std::string GetName() const { return "FiniteElement"; }

  int n() const { return B.n(); }
  double N(int i) const { return B.phi(i); }
  double N_x(int i) const { return grad[i].x(); }
  double N_y(int i) const { return grad[i].y(); }
  double N_z(int i) const { return grad[i].z(); }
  double J() const { return det; }
  double G() const { return T.G(); }

  void x(Vertex<DIM>& v) const { v = T.x(); }

  void mult_ad(Vertex<DIM>& p1, Vertex<DIM>& p2) const
  {
    T.DTI().mult_ad(p1, p2);
  }
  void normal(Vertex<DIM>& v) const { v = T.normal(); };

  void point(const Vertex<DIM>&) const;
  void point_T(const Vertex<DIM>& xi) const { T.point(xi); }
  void point_boundary(int ie, const Vertex<BDIM>& s1) const;
  /// depreciated
  void ReInit(const Matrix& M) const
  {
    assert(M.n() == DIM);
    assert(M.m() == B.n());
    T.ReInit(M);
  }

  void init_test_functions(TestFunction& Phi, double w, int i) const;

  void Anisotropy(DoubleMatrix& A) const;

  void GetCoordinates(DoubleMatrix& A) const { T.GetCoordinates(A); }
};

typedef FiniteElement<2, 1, Transformation2d<BaseQ12d>, BaseQ12d>
  FiniteElementQ12d;
typedef FiniteElement<2, 1, Transformation2d<BaseQ22d>, BaseQ22d>
  FiniteElementQ22d;
typedef FiniteElement<3, 2, Transformation3d<BaseQ13d>, BaseQ13d>
  FiniteElementQ13d;
typedef FiniteElement<3, 2, Transformation3d<BaseQ23d>, BaseQ23d>
  FiniteElementQ23d;

} // namespace Gascoigne

/*-----------------------------------------------------*/

#endif
