/*----------------------------   heatproblem.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __heatproblem_H
#define __heatproblem_H
/*----------------------------   heatproblem.h     ---------------------------*/

/**
 *
 * Copyright (C) 2020 by the Gascoigne 3D authors
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

#include "boundaryequation.h"
#include "dirichletdata.h"
#include "domainrighthandside.h"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"

/*-----------------------------------------*/

extern double __TIME;

namespace Gascoigne {

/**
 * Stores all the data required in the program that is
 * defined in the parameter file
 * This class is used in the equation but also within
 * Dirichlet data, right hand side, ...
 **/
class HeatData
{
public:
  double lambda, mu, rho, dt;

  void BasicInit(const ParamFile& pf)
  {
    DataFormatHandler DFH;
    DFH.insert("lambda", &lambda, 1.);
    DFH.insert("mu", &mu, 1.);
    DFH.insert("rho", &rho, 1.);
    DFH.insert("dt", &dt, 1.);
    FileScanner FS(DFH, pf, "Equation");
    assert(dt > 0);
    assert(mu > 0);
    assert(rho > 0);
    std::cout << dt << "\t" << mu << "\t" << lambda << "\t" << rho << std::endl;
  }
};

/**
 * Defines the DirichletData of the problem. T=0 everywhere
 **/
template<int DIM>
class HeatDD : public DirichletData
{
protected:
public:
  HeatDD(const ParamFile&)
    : DirichletData()
  {
    colors.clear();
    colors.insert(4);
    comp_on_color.clear();
    for (int i = 0; i < 2 * DIM; ++i)
      comp_on_color[4].push_back(i);
  }

  void operator()(DoubleVector& b,
                  const Vertex2d& /*unused*/,
                  int /*unused*/) const
  {
    b.zero();
  }
  void operator()(DoubleVector& b,
                  const Vertex3d& /*unused*/,
                  int /*unused*/) const
  {
    b.zero();
  }
};

/**
 * Defines the RHS of the problem. Depents on time (global variable)
 **/
template<int DIM>
class HeatRhs : public DomainRightHandSide
{
protected:
public:
  int GetNcomp() const { return 2 * DIM; }

  double operator()(int c, const Vertex2d& v) const
  {
    if (c != DIM - 1) // Kraft nur in y-Richtung (2d) und z-Richtung (3d)
      return 0.;

    double t = __TIME; // jeweils  Impuls fuer 1/10 sek.
    while (t > 1.0)
      t -= 1.0;
    if (t < 0.1)
      return -1.0;
    else
      return 0.0;
  }
  double operator()(int c, const Vertex3d& v) const
  {
    if (c != DIM - 1) // Kraft nur in y-Richtung (2d) und z-Richtung (3d)
      return 0.;

    double t = __TIME; // jeweils  Impuls fuer 1/10 sek.
    while (t > 1.0)
      t -= 1.0;
    if (t < 0.1)
      return -1.0;
    else
      return 0.0;
  }

  HeatRhs* createNew() const { return new HeatRhs<DIM>; }
};

template<int DIM>
class HeatEquation : public virtual Equation
{
  mutable FemFunction* old;

protected:
  HeatData data;

public:
  HeatEquation(const HeatData& PD)
    : data(PD)
  {}

  int GetNcomp() const { return 2 * DIM; } //

  void SetFemData(FemData& q) const
  {
    assert(q.find("OLD") != q.end());
    old = &q["OLD"];
  }

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    for (size_t i = 0; i < DIM; ++i) {
      // dt v
      b[i] += (U[i + DIM].m() - (*old)[i + DIM].m()) / data.dt * N.m();

      for (size_t j = 0; j < DIM; ++j) {
        b[i] +=
          2.0 * data.mu / data.rho * (U[i][j + 1] + U[j][i + 1]) * N[j + 1];
        b[i] += data.lambda / data.rho * U[j][j + 1] * N[i + 1];
      }

      // dt u - v = 0
      b[i + DIM] += (U[i].m() - (*old)[i].m()) / data.dt * N.m();
      b[i + DIM] -= U[i + DIM].m() * N.m();
    }
  }

  void Matrix(EntryMatrix& A,
              const FemFunction& /*unused*/,
              const TestFunction& M,
              const TestFunction& N) const
  {

    for (int i = 0; i < DIM; ++i) {
      // dt v
      A(i, i + DIM) += M.m() / data.dt * N.m();

      for (int j = 0; j < DIM; ++j) {
        A(i, i) += 2.0 * data.mu / data.rho * M[j + 1] * N[j + 1];
        A(i, j) += 2.0 * data.mu / data.rho * M[i + 1] * N[j + 1];
        A(i, j) += data.lambda / data.rho * M[j + 1] * N[i + 1];
      }

      // dt u - v = 0
      A(i + DIM, i) += M.m() / data.dt * N.m();
      A(i + DIM, i + DIM) -= M.m() * N.m();
    }
  }

  HeatEquation* createNew() const { return new HeatEquation<DIM>(data); }
};

} // namespace Gascoigne

/*----------------------------   heatproblem.h     ---------------------------*/
/* end of #ifndef __heatproblem_H */
#endif
/*----------------------------   heatproblem.h     ---------------------------*/
