#include "problemdescriptorbase.h"

namespace Gascoigne{

/*------------------------------------------------------------------------------*/

ProblemDescriptorBase::  ProblemDescriptorBase() : EQ(NULL),BM(NULL),ES(NULL),IC(NULL),
						   RHS(NULL),DD(NULL),BRHS(NULL),BE(NULL),
						   _paramfile(NULL) 
{}

/*------------------------------------------------------------------------------*/

ProblemDescriptorBase::~ProblemDescriptorBase() 
{
  if (EQ!=NULL)   { delete EQ;  EQ=NULL;}
  if (BM!=NULL)   { delete BM;  BM=NULL;}
  if (ES!=NULL)   { delete ES;  ES=NULL;}
  if (IC!=NULL)   { delete IC;  IC=NULL;}
  if (RHS!=NULL)  { delete RHS; RHS=NULL;}
  if (DD!=NULL)   { delete DD;  DD=NULL;}
  if (BRHS!=NULL) { delete BRHS;  BRHS=NULL;}
  if (BE!=NULL)   { delete BE;  BE=NULL;}
}

/*------------------------------------------------------------------------------*/

std::ostream& ProblemDescriptorBase::OutputSettings(std::ostream& os) const 
{
  if(EQ)  os << "Equation:        " << EQ->GetName() << std::endl;
  if(BM)  os << "BoundaryManager: " << BM->GetName() << std::endl;
  if(RHS) os << "Rhs:             " << RHS->GetName() << std::endl;
  if(DD)  os << "DirichletData:   " << DD->GetName() << std::endl;
  if(ES)  os << "ExactSolution:   " << ES->GetName() << std::endl;
  if(IC)  os << "InitialCondition:" << IC->GetName() << std::endl;
  return os;
}
  
/*------------------------------------------------------------------------------*/

void ProblemDescriptorBase::BasicInit(const ParamFile* pf) 
{
  if(GetParamFile()==NULL)
    {
      GetParamFilePointer() = pf;
    }
  if(GetBoundaryManagerPointer()==NULL)
    {
      GetBoundaryManagerPointer() = new BoundaryManager();
    }
  GetBoundaryManager()->BasicInit(_paramfile);
}

/*------------------------------------------------------------------------------*/

void ProblemDescriptorBase::SetTime(double time, double dt) const 
{
  if (EQ)  EQ  -> SetTime(time,dt);
  if (ES)  ES  -> SetTime(time,dt);
  if (RHS) RHS -> SetTime(time,dt);
  if (DD)  DD  -> SetTime(time,dt);
  if (BRHS) BRHS -> SetTime(time,dt);
  if (BE)   BE   -> SetTime(time,dt);
  if (IC)  IC  -> SetTime(time,dt);
}
}
