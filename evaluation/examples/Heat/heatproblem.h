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
  double lambda, dt;

  void BasicInit(const ParamFile& pf)
  {
    DataFormatHandler DFH;
    DFH.insert("lambda", &lambda, 1.);
    DFH.insert("dt", &dt, 1.);
    FileScanner FS(DFH, pf, "Equation");
    assert(lambda > 0);
    assert(dt > 0);
  }
};

/**
 * Defines the DirichletData of the problem. T=0 everywhere
 **/
class HeatDD : public DirichletData
{
protected:
public:
  HeatDD(const ParamFile& pf)
    : DirichletData(pf)
  {}

  void operator()(DoubleVector& b,
                  const Vertex2d& /*unused*/,
                  int /*unused*/) const
  {
    b.zero();
  }
};

/**
 * Defines the RHS of the problem. Depents on time (global variable)
 **/
class HeatRhs : public DomainRightHandSide
{
protected:
public:
  int GetNcomp() const { return 1.; }

  double operator()(int /*unused*/, const Vertex2d& v) const
  {
    double mx = 0.5 + 0.25 * cos(2.0 * M_PI * __TIME) - v.x();
    double my = 0.5 + 0.25 * sin(2.0 * M_PI * __TIME) - v.y();
    return exp(-25.0 * (mx * mx + my * my));
  }

  HeatRhs* createNew() const { return new HeatRhs; }
};

class HeatEquation : public virtual Equation
{
  mutable FemFunction* old;

protected:
  HeatData data;

public:
  HeatEquation(const HeatData& PD)
    : data(PD)
  {}

  int GetNcomp() const { return 1; } // only one component

  void SetFemData(FemData& q) const
  {
    assert(q.find("OLD") != q.end());
    old = &q["OLD"];
  }

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    b[0] += (U[0].m() - (*old)[0].m()) / data.dt * N.m();
    b[0] += data.lambda * (U[0].x() * N.x() + U[0].y() * N.y());
    b[0] += -U[0].y() * N.m();
  }

  void Matrix(EntryMatrix& A,
              const FemFunction& /*unused*/,
              const TestFunction& M,
              const TestFunction& N) const
  {
    A(0, 0) += M.m() / data.dt * N.m();
    A(0, 0) += data.lambda * (M.x() * N.x() + M.y() * N.y());
    A(0, 0) += -M.y() * N.m();
  }

  HeatEquation* createNew() const { return new HeatEquation(data); }
};

} // namespace Gascoigne

/*----------------------------   heatproblem.h     ---------------------------*/
/* end of #ifndef __heatproblem_H */
#endif
/*----------------------------   heatproblem.h     ---------------------------*/
