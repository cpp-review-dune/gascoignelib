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

#include "problemdescriptorbase.h"
#include "componentinformationbase.h"

using namespace std;

namespace Gascoigne {

/*------------------------------------------------------------------------------*/

ProblemDescriptorBase::ProblemDescriptorBase()
  : EQ(NULL)
  , BEQ(NULL)
  , RHS(NULL)
  , DRHS(NULL)
  , BRHS(NULL)
  , FEQ(NULL)
  , BM(NULL)
  , ES(NULL)
  , IC(NULL)
  , DD(NULL)
  , PD(NULL)
  , BIC(NULL)
  , CI(NULL)
  , _time(0.0)
  , _dt(0.0)
{}

/*------------------------------------------------------------------------------*/

ProblemDescriptorBase::~ProblemDescriptorBase()
{
  if (BM != NULL) {
    delete BM;
    BM = NULL;
  }
  if (ES != NULL) {
    delete ES;
    ES = NULL;
  }
  if (IC != NULL) {
    delete IC;
    IC = NULL;
  }
  if (DD != NULL) {
    delete DD;
    DD = NULL;
  }
  if (PD != NULL) {
    delete PD;
    PD = NULL;
  }
  if (BIC != NULL) {
    delete BIC;
    BIC = NULL;
  }
  if (CI != NULL) {
    delete CI;
    CI = NULL;
  }
}

/*------------------------------------------------------------------------------*/

ostream&
ProblemDescriptorBase::OutputSettings(ostream& os) const
{
  if (DD)
    os << "DirichletData:            " << DD->GetName() << endl;
  if (PD)
    os << "PeriodicData:             " << PD->GetName() << endl;
  if (ES)
    os << "ExactSolution:            " << ES->GetName() << endl;
  if (IC)
    os << "InitialCondition:         " << IC->GetName() << endl;
  if (BIC)
    os << "BoundaryInitialCondition: " << BIC->GetName() << endl;
  if (BM)
    os << "BoundaryManager:          " << BM->GetName() << endl;
  if (CI)
    os << "ComponentInformation:     " << CI->GetName() << std::endl;
  return os;
}

/*------------------------------------------------------------------------------*/

void
ProblemDescriptorBase::BasicInit(const ParamFile& pf)
{
  _paramfile = pf;

  // read all input to control the solver
  GetSolverData().BasicInit(pf);

  if (GetBoundaryManagerPointer() == NULL) {
    GetBoundaryManagerPointer() = new BoundaryManager();
  }
  GetBoundaryManager()->BasicInit(_paramfile);

  if (GetComponentInformation() == NULL) {
    GetComponentInformationPointer() = new ComponentInformationBase();
  }
  GetComponentInformationPointer()->GetProblemDescriptorInterface() = this;
}

/*------------------------------------------------------------------------------*/

void
ProblemDescriptorBase::SetTime(double time, double dt) const
{
  _time = time;
  _dt = dt;

  if (FEQ)
    FEQ->SetTime(time, dt);
  if (BRHS)
    BRHS->SetTime(time, dt);
  if (DRHS)
    DRHS->SetTime(time, dt);
  if (RHS)
    RHS->SetTime(time, dt);
  if (BEQ)
    BEQ->SetTime(time, dt);
  if (EQ)
    EQ->SetTime(time, dt);
  if (ES)
    ES->SetTime(time, dt);
  if (DD)
    DD->SetTime(time, dt);
  if (PD)
    PD->SetTime(time, dt);
  if (IC)
    IC->SetTime(time, dt);
  if (BIC)
    BIC->SetTime(time, dt);
  if (CI)
    CI->SetTime(time, dt);
}
} // namespace Gascoigne
