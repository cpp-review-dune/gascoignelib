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
#include <Solver/stdloop.h>

#include "heatproblem.h"

extern double __TIME;

namespace Gascoigne {

template<int COMP>
class HeatProblem : public Gascoigne::ProblemDescriptorBase
{
  Gascoigne::HeatData data;

public:
  void BasicInit(const Gascoigne::ParamFile& pf)
  {
    data.BasicInit(pf);

    GetParamFile() = pf;
    GetEquationPointer() = new Gascoigne::HeatEquation<COMP>(data);
    GetDirichletDataPointer() = new Gascoigne::HeatDD<COMP>(pf);
    GetRightHandSidePointer() = new Gascoigne::HeatRhs<COMP>;

    ProblemDescriptorBase::BasicInit(pf);
  }
};
} // namespace Gascoigne

#endif
