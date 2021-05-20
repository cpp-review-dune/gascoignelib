/**
 *
 * Copyright (C) 2004, 2007 by the Gascoigne 3D authors
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

#ifndef __LOCAL1_h
#define __LOCAL1_h

#include "dgequation.h"
#include "domainmeanfunctional.h"
#include "equation.h"
#include "problemdescriptorbase.h"
#include "residualfunctional.h"
#include "righthandsidebyequation.h"
#include "zerodirichletdata.h"

/*---------------------------------------------------*/

namespace Gascoigne {

class MyEQ : public DGEquation
{
  mutable double f;

  mutable bool _internaledge;
  mutable double _h;
  mutable Vertex2d _n;

public:
  int GetNcomp() const { return 1; }
  std::string GetName() const { return "DGEQ"; }
  void point(double h, const FemFunction& U, const Vertex2d& v) const
  {
    f = 2.0 * M_PI * M_PI * sin(M_PI * v.x()) * sin(M_PI * v.y());
  }

  void point_edge(bool internaledge,
                  double h,
                  const FemFunction& U1,
                  const FemFunction& U2,
                  const Vertex2d& v,
                  const Vertex2d& n) const
  {
    _internaledge = internaledge;
    _h = h;
    _n = n;
  }
  void EdgeForm1(VectorIterator b,
                 const FemFunction& U1,
                 const FemFunction& U2,
                 const TestFunction& N) const
  {
    // continuity and boundary
    if (_internaledge)
      b[0] += -1.e2 / _h * (U2[0].m() - U1[0].m()) * N.m();
    else
      b[0] += -1.e2 / _h * (-U1[0].m()) * N.m();

    // repair normal
    if (_internaledge) {
      b[0] -= 0.5 * (_n.x() * U1[0].x() + _n.y() * U1[0].y()) * N.m();
      b[0] -= 0.5 * (_n.x() * U2[0].x() + _n.y() * U2[0].y()) * N.m();
    } else {
      b[0] -= (_n.x() * U1[0].x() + _n.y() * U1[0].y()) * N.m();
    }
  }
  void EdgeForm2(VectorIterator b,
                 const FemFunction& U1,
                 const FemFunction& U2,
                 const TestFunction& N) const
  {
    // continuity
    b[0] += 1.e2 / _h * (U2[0].m() - U1[0].m()) * N.m();

    // repair normal - other side, reverse normal? Jump?
    b[0] += 0.5 * (_n.x() * U1[0].x() + _n.y() * U1[0].y()) * N.m();
    b[0] += 0.5 * (_n.x() * U2[0].x() + _n.y() * U2[0].y()) * N.m();
  }
  void EdgeMatrix11(EntryMatrix& A,
                    const FemFunction& U1,
                    const FemFunction& U2,
                    const TestFunction& M,
                    const TestFunction& N) const
  {
    // continuity and boundary
    A(0, 0) += -1.e2 / _h * (-M.m()) * N.m();

    // repair normal
    if (_internaledge) {
      A(0, 0) -= 0.5 * (_n.x() * M.x() + _n.y() * M.y()) * N.m();
    } else {
      A(0, 0) -= (_n.x() * M.x() + _n.y() * M.y()) * N.m();
    }
  }

  void EdgeMatrix12(EntryMatrix& A,
                    const FemFunction& U1,
                    const FemFunction& U2,
                    const TestFunction& M,
                    const TestFunction& N) const
  {
    // continuity and boundary
    A(0, 0) += -1.e2 / _h * (M.m()) * N.m();

    // repair normal
    A(0, 0) -= 0.5 * (_n.x() * M.x() + _n.y() * M.y()) * N.m();
  }
  void EdgeMatrix21(EntryMatrix& A,
                    const FemFunction& U1,
                    const FemFunction& U2,
                    const TestFunction& M,
                    const TestFunction& N) const
  {
    // continuity and boundary
    A(0, 0) += 1.e2 / _h * (-M.m()) * N.m();

    // repair normal
    A(0, 0) += 0.5 * (_n.x() * M.x() + _n.y() * M.y()) * N.m();
  }
  void EdgeMatrix22(EntryMatrix& A,
                    const FemFunction& U1,
                    const FemFunction& U2,
                    const TestFunction& M,
                    const TestFunction& N) const
  {
    // continuity and boundary
    A(0, 0) += 1.e2 / _h * (M.m()) * N.m();

    // repair normal
    A(0, 0) += 0.5 * (_n.x() * M.x() + _n.y() * M.y()) * N.m();
  }

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    b[0] += U[0].x() * N.x() + U[0].y() * N.y();
    b[0] += -f * N.m();
  }
  void Matrix(EntryMatrix& A,
              const FemFunction& U,
              const TestFunction& M,
              const TestFunction& N) const
  {
    A(0, 0) += M.x() * N.x() + M.y() * N.y();
  }
};

class ProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
public:
  std::string GetName() const { return "Local"; }
  void BasicInit(const Gascoigne::ParamFile* pf)
  {
    GetEquationPointer() = new MyEQ;
    //      GetRightHandSidePointer() = new MyRhs;
    GetDirichletDataPointer() = new ZeroDirichletData();

    ProblemDescriptorBase::BasicInit(pf);
  }
};
} // namespace Gascoigne

#endif
