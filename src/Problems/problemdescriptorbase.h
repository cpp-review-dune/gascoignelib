#ifndef __ProblemDescriptorBase_h
#define __ProblemDescriptorBase_h

#include "problemdescriptorinterface.h"

/**********************************************************/

namespace Gascoigne
{
class ProblemDescriptorBase : public ProblemDescriptorInterface
{

 private:
  
  Equation              *EQ;
  BoundaryManager       *BM;
  ExactSolution         *ES;
  Application           *RHS;
  Application           *IC;
  DirichletData         *DD;
  BoundaryRightHandSide *BRHS;
  BoundaryEquation      *BE;
  
  const ParamFile *_paramfile;
  
 protected:
  
  const ParamFile*& GetParamFilePointer() {return _paramfile;}
  
  Equation*& GetEquationPointer() { return EQ;}
  BoundaryManager*& GetBoundaryManagerPointer() { return BM;}
  ExactSolution*& GetExactSolutionPointer() { return ES;}
  Application*& GetRightHandSidePointer() { return RHS;}
  Application*& GetInitialConditionPointer() { return IC;}
  DirichletData*& GetDirichletDataPointer() { return DD;}
  BoundaryRightHandSide*& GetBoundaryRightHandSidePointer() { return BRHS;}
  BoundaryEquation*& GetBoundaryEquationPointer() { return BE; }

  BoundaryManager* GetBoundaryManager() { return BM;}
  
 public:

  ProblemDescriptorBase() : EQ(NULL),BM(NULL),ES(NULL),IC(NULL),RHS(NULL),DD(NULL),BRHS(NULL),BE(NULL),_paramfile(NULL) {}
  ~ProblemDescriptorBase() {
    if (EQ!=NULL)   { delete EQ;  EQ=NULL;}
    if (BM!=NULL)   { delete BM;  BM=NULL;}
    if (ES!=NULL)   { delete ES;  ES=NULL;}
    if (IC!=NULL)   { delete IC;  IC=NULL;}
    if (RHS!=NULL)  { delete RHS; RHS=NULL;}
    if (DD!=NULL)   { delete DD;  DD=NULL;}
    if (BRHS!=NULL) { delete BRHS;  BRHS=NULL;}
    if (BE!=NULL)   { delete BE;  BE=NULL;}
  }

  std::ostream& OutputSettings(std::ostream& os) const {
    if(EQ)   os << "Equation:        " << EQ->GetName()   << std::endl;
    if(BE)   os << "BoundaryEquation:" << BE->GetName()   << std::endl;
    if(RHS)  os << "Rhs:             " << RHS->GetName()  << std::endl;
    if(BRHS) os << "BoundaryRhs:     " << BRHS->GetName() << std::endl;
    if(DD)   os << "DirichletData:   " << DD->GetName()   << std::endl;
    if(ES)   os << "ExactSolution:   " << ES->GetName()   << std::endl;
    if(IC)   os << "InitialCondition:" << IC->GetName()   << std::endl;
    if(BM)   os << "BoundaryManager: " << BM->GetName()   << std::endl;
    return os;
  }
  
  void BasicInit(const ParamFile* pf) {
    if(!GetParamFile())
      {
	GetParamFilePointer() = pf;
      }
    if(!GetBoundaryManagerPointer())
      {
	GetBoundaryManagerPointer() = new BoundaryManager();
      }
    GetBoundaryManager()->BasicInit(_paramfile);
  }
  
  const ParamFile* GetParamFile() const {return _paramfile;}

  const Application*           GetRightHandSide        () const { return  RHS;}
  const DirichletData*         GetDirichletData        () const { return  DD;}
  const BoundaryRightHandSide* GetBoundaryRightHandSide() const { return  BRHS;}
  const BoundaryEquation*      GetBoundaryEquation     () const { return  BE;}
  const Application*           GetInitialCondition     () const { return IC;}
  const ExactSolution*         GetExactSolution        () const { return ES;}
  const Equation*              GetEquation             () const { return EQ;}
  const BoundaryManager*       GetBoundaryManager      () const { return BM;}

  void SetTime(double time, double dt) const {
    if (EQ)   EQ   -> SetTime(time,dt);
    if (ES)   ES   -> SetTime(time,dt);
    if (RHS)  RHS  -> SetTime(time,dt);
    if (DD)   DD   -> SetTime(time,dt);
    if (BRHS) BRHS -> SetTime(time,dt);
    if (BE)   BE   -> SetTime(time,dt);
    if (IC)   IC   -> SetTime(time,dt);
  }
};
}

/**********************************************************/

#endif
