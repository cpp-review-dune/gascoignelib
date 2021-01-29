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


#ifndef __ProblemDescriptorBase_h
#define __ProblemDescriptorBase_h

#include "problemdescriptorinterface.h"

/*--------------------------------------------------*/

namespace Gascoigne
{
  class ProblemDescriptorBase : public ProblemDescriptorInterface
  {
  public:


    
  private:
    Equation* EQ;
    BoundaryEquation* BEQ;
    DomainRightHandSide* RHS;
    DiracRightHandSide* DRHS;
    BoundaryRightHandSide* BRHS;

    FaceEquation *FEQ;
    BoundaryManager *BM;
    ExactSolution *ES;
    Application *IC;
    DirichletData *DD;
    PeriodicData *PD;
    BoundaryInitialCondition *BIC;
    ComponentInformation *CI;

    ParamFile _paramfile;
    mutable double _time,_dt;

  protected:

    Equation*& GetEquationPointer()
    {
      return EQ;
    }
    BoundaryEquation*& GetBoundaryEquationPointer()
    {
      return BEQ;
    }
    DomainRightHandSide*& GetRightHandSidePointer()
    {
      return RHS;
    }
    DiracRightHandSide*& GetDiracRightHandSidePointer()
    {
      return DRHS;
    }
    BoundaryRightHandSide*& GetBoundaryRightHandSidePointer()
    {
      return BRHS;
    }
    
    FaceEquation *&GetFaceEquationPointer()
    {
      return FEQ;
    }
    BoundaryManager *&GetBoundaryManagerPointer()
    {
      return BM;
    }
    ExactSolution *&GetExactSolutionPointer()
    {
      return ES;
    }
    Application *&GetInitialConditionPointer()
    {
      return IC;
    }
    DirichletData *&GetDirichletDataPointer()
    {
      return DD;
    }
    PeriodicData *&GetPeriodicDataPointer()
    {
      return PD;
    }
    BoundaryInitialCondition *&GetBoundaryInitialConditionPointer()
    {
      return BIC;
    }
    BoundaryManager *GetBoundaryManager()
    {
      return BM;
    }

    ComponentInformation *&GetComponentInformationPointer()
    {
      return CI;
    }

  public:
    ProblemDescriptorBase();

    ~ProblemDescriptorBase();

    std::ostream &OutputSettings(std::ostream &os) const;

    void BasicInit(const ParamFile &pf);

  // Gives the number of solution components. 
    int GetNcomp() const 
    {
      return GetEquation()->GetNcomp();
    }
    

    const ParamFile& GetParamFile() const { return _paramfile; }
    ParamFile& GetParamFile() { return _paramfile; }

    const Equation* GetEquation() const
    {
      return EQ;
    }
    const BoundaryEquation* GetBoundaryEquation() const
    {
      return BEQ;
    }
    const DomainRightHandSide* GetRightHandSide() const
    {
      return RHS;
    }
    const DiracRightHandSide* GetDiracRightHandSide() const
    {
      return DRHS;
    }
    const BoundaryRightHandSide* GetBoundaryRightHandSide() const
    {
      return BRHS;
    }
    
    const DirichletData *GetDirichletData() const
    {
      return DD;
    }
    const PeriodicData *GetPeriodicData() const
    {
      return PD;
    }
    const BoundaryInitialCondition *GetBoundaryInitialCondition() const
    {
      return BIC;
    }
    const Application *GetInitialCondition() const
    {
      return IC;
    }
    const ExactSolution *GetExactSolution() const
    {
      return ES;
    }
    const FaceEquation *GetFaceEquation() const
    {
      return FEQ;
    }
    const BoundaryManager *GetBoundaryManager() const
    {
      return BM;
    }
    const ComponentInformation *GetComponentInformation() const
    {
      return CI;
    }

    // stores the time and the time step
    // (should be moved to the problem data)
    void SetTime(double time, double dt) const;
    // access to time
    double time() const { return _time; }
    // access to time step
    double dt() const   { return _dt; }     
    
  };


} // namespace Gascoigne

/*--------------------------------------------------*/

#endif
