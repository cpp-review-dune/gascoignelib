/**
 *
 * Copyright (C) 2004, 2005, 2007 by the Gascoigne 3D authors
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

#ifndef __local_h
#define __local_h

#include <Common/stopwatch.h>
#include <Problems/problemdescriptorbase.h>
#include <Solver/componentinformationbase.h>
#include <Solver/stdloop.h>

#include "elasticityproblem.h"

extern double __TIME;

namespace Gascoigne {

template<IndexType DIM>
class MyCI : public ComponentInformationBase
{
public:
  const IndexType GetNScalars() const { return 2 * DIM; }

  void GetScalarName(IndexType i, std::string& s_name) const
  {
    if (DIM == 2) {
      if (i == 0)
        s_name = "Ux";
      else if (i == 1)
        s_name = "Uy";
      else if (i == 2)
        s_name = "Vx";
      else if (i == 3)
        s_name = "Vy";
      else
        abort();
    } else if (DIM == 3) {
      if (i == 0)
        s_name = "Ux";
      else if (i == 1)
        s_name = "Uy";
      else if (i == 2)
        s_name = "Uz";
      else if (i == 3)
        s_name = "Vx";
      else if (i == 4)
        s_name = "Vy";
      else if (i == 5)
        s_name = "Vz";
      else
        abort();
    } else
      abort();
  }

  virtual const IndexType GetNVectors() const { return 2; }

  virtual void GetVectorName(IndexType i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "U";
    else
      s_name = "V";
  }

  virtual void GetVectorIndices(
    IndexType i,
    std::array<IndexType, 3>& fa_vectorindices) const
  {
    if (i == 0)
      fa_vectorindices = { 0, 1, 2 };
    else
      fa_vectorindices = { DIM, DIM + 1, DIM + 2 };
    if (DIM == 2)
      fa_vectorindices[2] = -1;
  }
};

template<int DIM>
class HeatProblem : public Gascoigne::ProblemDescriptorBase
{
  Gascoigne::HeatData data;

public:
  Gascoigne::HeatData& GetData() { return data; }

  void BasicInit(const Gascoigne::ParamFile& pf)
  {
    data.BasicInit(pf);

    GetParamFile() = pf;
    GetEquationPointer() = new Gascoigne::HeatEquation<DIM>(data);
    GetDirichletDataPointer() = new Gascoigne::HeatDD<DIM>(pf);
    GetRightHandSidePointer() = new Gascoigne::HeatRhs<DIM>;
    GetComponentInformationPointer() = new Gascoigne::MyCI<DIM>;

    ProblemDescriptorBase::BasicInit(pf);
  }
};
} // namespace Gascoigne

#endif
