/**
 *
 * Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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

#include "gmres.h"
#include "compose_name.h"
#include "nmatrix.h"

namespace Gascoigne {

/********************************************************************/

GMRES::GMRES(StdSolver& S_, StdMultiLevelSolver& P_, int vm)
  : H(vm)
  , gamma(vm + 1)
  , ci(vm)
  , si(vm)
  , vmax(vm)
  , S(S_)
  , P(P_)
{
  left_precondition = 1;
  for (int i = 0; i < 3; i++) {
    new_memory();
  }
}

/********************************************************************/

GMRES::~GMRES()
{
  for (int i = 0; i < mem.size(); ++i) {
    P.DeleteVector(mem[i]);
  }
  mem.clear();
}

/********************************************************************/

void
GMRES::new_memory()
{
  int i = mem.size();

  std::string s = "gmres";
  compose_name(s, i);
  mem.resize(i + 1, s);
  // mem[i].SetMultiLevelSolver(&P);
  P.ReInitVector(mem[i]);
}

/********************************************************************/

void
GMRES::init()
{
  assert(H.size() == vmax);
  for (int i = 0; i < H.size(); i++) {
    H[i].zero();
    H[i].resize(i + 2, 0.);
  }
  gamma.zero();
  ci.zero();
  si.zero();
}

/********************************************************************/

void
GMRES::givens_rotation(dvector& h, int col)
{
  /*  Transformation into triagonal structure  */

  for (int i = 0; i < col; i++) {
    double dummy = h[i];
    h[i] = ci[i] * dummy + si[i] * h[i + 1];
    h[i + 1] = -si[i] * dummy + ci[i] * h[i + 1];
  }
  double r = 1. / sqrt(h[col] * h[col] + h[col + 1] * h[col + 1]);
  si[col] = h[col + 1] * r;
  ci[col] = h[col] * r;
  h[col] = ci[col] * h[col] + si[col] * h[col + 1];
  gamma[col + 1] = -si[col] * gamma[col];
  gamma[col] *= ci[col];
}

/********************************************************************/

double
GMRES::orthogonalization(dvector& h, int dim, Vector& vv) const
{
  for (int i = 0; i < dim; i++) {
    double d = S.ScalarProduct(vv, mem[i]);
    h[i] += d;
    S.Add(vv, -d, mem[i]);
  }
  h[dim] = sqrt(S.ScalarProduct(vv, vv));
  return h[dim];
}

/********************************************************************/

bool
GMRES::reortho_test(const Vector& u, double norm) const
{
  bool test = 0;
  double delta = 1.e-3;
  double oldnorm = sqrt(S.ScalarProduct(u, u));
  double newnorm = oldnorm + delta * norm;

  if (oldnorm == newnorm) {
    std::cerr << "REORTHO !!\n";
    test = 1;
  }
  return test;
}

/********************************************************************/

int
GMRES::restarted_solve(const Matrix& A,
                       Vector& x,
                       const Vector& b,
                       CGInfo& info)
{
  int reached;
  for (int j = 0;; j++) {
    info.control().iteration() = j * (vmax - 1);
    reached = solve(A, x, b, info);
    if (reached)
      break;
  }
  if (reached < 0)
    return 1;
  return 0;
}

/********************************************************************/

int
GMRES::solve(const Matrix& A, Vector& x, const Vector& b, CGInfo& info)
{
  init();

  int k0 = info.control().iteration();

  int reached = 0;

  Vector& v = mem[0];

  Vector p("gmresp");
  // p.SetMultiLevelSolver(&P);
  P.ReInitVector(p);
  P.Equ(p, 1., v);

  if (left_precondition) {
    S.residualgmres(A, p, x, b);
    P.precondition(A, v, p);
  } else {
    S.residualgmres(A, v, x, b);
  }

  double rho = sqrt(S.ScalarProduct(v, v));
  if (rho == 0.) {
    P.DeleteVector(p);
    return 1;
  }

  gamma[0] = rho;

  // reached = info.check(rho,0);

  S.Equ(v, 1. / rho, v);

  int dim = 0;
  while ((dim < vmax - 1) && !reached) {
    dim++;
    int last = dim - 1;
    if (dim >= mem.size()) {
      new_memory();
    }
    Vector& vv = mem[dim];
    if (left_precondition) {
      S.vmulteq(A, p, mem[last], 1.);
      P.precondition(A, vv, p);
    } else {
      S.Zero(p);
      P.precondition(A, p, mem[last]);
      S.vmulteq(A, vv, p, 1.);
    }
    double s = orthogonalization(H[last], dim, vv);
    if (reortho_test(vv, s)) {
      s = orthogonalization(H[last], dim, vv);
    }
    S.Equ(vv, 1. / s, vv);

    givens_rotation(H[last], last);

    /*  residual  */

    rho = fabs(gamma[dim]);

    reached = info.check(rho, k0 + last);
  }
  /*  Calculate solution  */

  solution(A, x, p, dim);

  P.DeleteVector(p);
  return reached;
}

/********************************************************************/

void
GMRES::solution(const Matrix& A, Vector& x, Vector& p, int dim)
{
  nmatrix<double> H1(dim + 1, dim);
  H1.zero();
  for (int j = 0; j < dim; j++) {
    for (int i = 0; i < j + 2; i++)
      H1(i, j) = H[j][i];
  }

  nvector<double> h(dim);

  backward(h, H1, gamma);

  if (left_precondition) {
    for (int i = 0; i < dim; i++)
      S.Add(x, h[i], mem[i]);
  } else {
    S.Zero(p);
    for (int i = 0; i < dim; i++)
      S.Add(p, h[i], mem[i]);
    S.Equ(mem[0], 0., mem[0]);
    P.precondition(A, mem[0], p);
    S.Add(x, 1., mem[0]);
  }
}
} // namespace Gascoigne
