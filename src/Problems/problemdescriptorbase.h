#ifndef __ProblemDescriptorBase_h
#define __ProblemDescriptorBase_h

#include "problemdescriptorinterface.h"

/**********************************************************/

namespace Gascoigne
{
class ProblemDescriptorBase : public ProblemDescriptorInterface
{

 private:
  
  Equation           *EQ;
  BoundaryManager    *BM;
  ExactSolution      *ES;
  InitialCondition   *IC;
  RightHandSideData  *RHS;
  DirichletData      *DD;
  NeumannData        *ND;
  
  const ParamFile* _paramfile;
  
 protected:
  
  const ParamFile* GetParamFile() const {return _paramfile;}
  
  Equation*& GetEquationPointer() { return EQ;}
  BoundaryManager*& GetBoundaryManagerPointer() { return BM;}
  ExactSolution*& GetExactSolutionPointer() { return ES;}
  InitialCondition*& GetInitialConditionPointer() { return IC;}
  RightHandSideData*& GetRightHandSideDataPointer() { return RHS;}
  DirichletData*& GetDirichletDataPointer() { return DD;}
  NeumannData*& GetNeumannDataPointer() { return ND;}

  BoundaryManager* GetBoundaryManager () { return BM;}

 public:

  ProblemDescriptorBase() : EQ(NULL),BM(NULL),ES(NULL),IC(NULL),RHS(NULL),DD(NULL),ND(NULL),_paramfile(NULL) {}
  ~ProblemDescriptorBase() {
    if (EQ!=NULL)   { delete EQ; EQ=NULL;}
    if (BM!=NULL)   { delete BM; BM=NULL;}
    if (ES!=NULL)   { delete ES; ES=NULL;}
    if (IC!=NULL)   { delete IC; IC=NULL;}
    if (RHS!=NULL)  { delete RHS; RHS=NULL;}
    if (DD!=NULL)   { delete DD; DD=NULL;}
    if (ND!=NULL)   { delete ND; ND=NULL;}
  }

  std::ostream& OutputSettings(std::ostream& os) const {
    if(EQ)  os << "Equation:        " << EQ->GetName() << std::endl;
    if(BM)  os << "BoundaryManager: " << BM->GetName() << std::endl;
    if(RHS) os << "Rhs:             " << RHS->GetName() << std::endl;
    if(DD)  os << "DirichletData:   " << DD->GetName() << std::endl;
    if(ND)  os << "NeumannData:     " << ND->GetName() << std::endl;
    if(ES)  os << "ExactSolution:   " << ES->GetName() << std::endl;
    if(IC)  os << "InitialCondition:" << IC->GetName() << std::endl;
    return os;
  }
  
  void BasicInit(const ParamFile* pf) {
    _paramfile = pf;
    if(GetBoundaryManagerPointer()==NULL)
      {
	GetBoundaryManagerPointer() = new BoundaryManager(_paramfile);
      }
  }
  
  const RightHandSideData*  GetRightHandSideData() const { return  RHS;}
  const DirichletData*      GetDirichletData()     const { return  DD;}
  const NeumannData*        GetNeumannData()       const { return  ND;}
  const InitialCondition*   GetInitialCondition() const { return IC;}
  const ExactSolution*   GetExactSolution   () const { return ES;}
  const Equation*        GetEquation        () const { return EQ;}
  const BoundaryManager* GetBoundaryManager () const { return BM;}

  void SetTime(double time, double dt) const {
    if (EQ)  EQ  -> SetTime(time,dt);
    if (ES)  ES  -> SetTime(time,dt);
    if (RHS) RHS -> SetTime(time,dt);
    if (DD)  DD  -> SetTime(time,dt);
    if (ND)  ND  -> SetTime(time,dt);
    if (IC)  IC  -> SetTime(time,dt);
  }
};
}

/**********************************************************/

#endif
