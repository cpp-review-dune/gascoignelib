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

#include "heatproblem.h"
#include "problemdescriptorbase.h"
#include "stdloop.h"
#include "stopwatch.h"

extern double __TIME;

namespace Gascoigne {

class HeatProblem : public Gascoigne::ProblemDescriptorBase
{
  Gascoigne::HeatData data;

public:
  void BasicInit(const Gascoigne::ParamFile& pf)
  {
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
    data.BasicInit(pf);

    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
    GetParamFile() = pf;
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
    GetEquationPointer() = new Gascoigne::HeatEquation(data);
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
    GetDirichletDataPointer() = new Gascoigne::HeatDD(pf);
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
    GetRightHandSidePointer() = new Gascoigne::HeatRhs;
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \

    ProblemDescriptorBase::BasicInit(pf);
  }
};
} // namespace Gascoigne

#endif
