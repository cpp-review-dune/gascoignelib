/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2009 by the Gascoigne 3D authors
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

#ifndef __ProblemDescriptorInterface_h
#define __ProblemDescriptorInterface_h

#include "filescanner.h"

#include "boundaryequation.h"
#include "boundarymanager.h"
#include "boundaryrighthandside.h"
#include "componentinformation.h"
#include "diracrighthandside.h"
#include "dirichletdata.h"
#include "domainrighthandside.h"
#include "equation.h"
#include "exactsolution.h"
#include "faceequation.h"
#include "periodicdata.h"

#include "solverdata.h"

namespace Gascoigne {
/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptorInterface

///
///
/////////////////////////////////////////////

class ProblemDescriptorInterface
{
private:
protected:
  SolverData solverdata;

public:
  ProblemDescriptorInterface() {}
  virtual ~ProblemDescriptorInterface() {}

  virtual void BasicInit(const ParamFile& pf) {}

  virtual std::string GetName() const { return "No Name"; }

  // Gives the number of solution components.
  virtual ShortIndexType GetNcomp() const
  {
    std::cerr << "ProblemDescriptor::GetNcomp() must be written!" << std::endl;
    abort();
  }

  virtual std::ostream& OutputSettings(std::ostream& os) const = 0;

  // stores the time and the time step
  // (should be moved to the problem data)
  virtual void SetTime(double time, double dt) const = 0;
  // access to time
  virtual double time() const
  {
    std::cerr << "PDI::time()" << std::endl;
    abort();
  }
  // access to time step
  virtual double dt() const
  {
    std::cerr << "PDI::time()" << std::endl;
    abort();
  }

  virtual const ParamFile& GetParamFile() const = 0;

  //// Description of the Numerics
  virtual const SolverData& GetSolverData() const { return solverdata; }
  virtual SolverData& GetSolverData() { return solverdata; }

  //// Description of the Problem
  virtual const Equation* GetEquation() const = 0;
  virtual const BoundaryEquation* GetBoundaryEquation() const = 0;
  virtual const DomainRightHandSide* GetRightHandSide() const = 0;
  virtual const DiracRightHandSide* GetDiracRightHandSide() const = 0;
  virtual const BoundaryRightHandSide* GetBoundaryRightHandSide() const = 0;
  virtual const FaceEquation* GetFaceEquation() const = 0;
  virtual const DirichletData* GetDirichletData() const = 0;
  virtual const PeriodicData* GetPeriodicData() const = 0;
  virtual const Application* GetInitialCondition() const = 0;
  virtual const BoundaryInitialCondition* GetBoundaryInitialCondition()
    const = 0;
  virtual const ExactSolution* GetExactSolution() const = 0;
  virtual const BoundaryManager* GetBoundaryManager() const = 0;
  virtual const ComponentInformation* GetComponentInformation() const = 0;
};
} // namespace Gascoigne

#endif
