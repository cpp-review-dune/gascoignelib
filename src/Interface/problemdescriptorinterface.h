#ifndef  __ProblemDescriptorInterface_h
#define  __ProblemDescriptorInterface_h

#include  "gostream.h"
#include  "stringutil.h"
#include  "filescanner.h"

#include  "boundarymanager.h"
#include  "equation.h"
#include  "righthandsidedata.h"
#include  "dirichletdata.h"
#include  "neumanndata.h"
#include  "exactsolution.h"
#include  "initialcondition.h"
#include  "boundarymanager.h"
#include  "data.h"

/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptorInterface

///
///
/////////////////////////////////////////////



class ProblemDescriptorInterface
{
public:


private:

  Equation           *EQ;
  BoundaryManager    *BM;
  ExactSolution      *ES;
  InitialCondition   *IC;
  RightHandSideData  *RHS;
  DirichletData      *DD;
  NeumannData        *ND;

  mutable std::string _paramfile;

protected:

  std::string GetParamFile() const {return _paramfile;}

  Equation*& GetEquationPointer() { return EQ;}
  BoundaryManager*& GetBoundaryManagerPointer() { return BM;}
  ExactSolution*& GetExactSolutionPointer() { return ES;}
  InitialCondition*& GetInitialConditionPointer() { return IC;}
  RightHandSideData*& GetRightHandSideDataPointer() { return RHS;}
  DirichletData*& GetDirichletDataPointer() { return DD;}
  NeumannData*& GetNeumannDataPointer() { return ND;}

  virtual void ConstructEquation()=0;
  virtual void ConstructRightHandSideData() {}
  virtual void ConstructDirichletData() {}
  virtual void ConstructNeumannData() {}
  virtual void ConstructExactSolution() {}
  virtual void ConstructInitialCondition() {}
  virtual void ConstructBoundaryManager() {
    BM = new BoundaryManager(_paramfile);
  }

public:

  ProblemDescriptorInterface() : EQ(NULL),BM(NULL),ES(NULL),IC(NULL),RHS(NULL),DD(NULL),ND(NULL)
{}
  virtual ~ProblemDescriptorInterface() {
    if (EQ!=NULL)   { delete EQ; EQ=NULL;}
    if (BM!=NULL)   { delete BM; BM=NULL;}
    if (ES!=NULL)   { delete ES; ES=NULL;}
    if (IC!=NULL)   { delete IC; IC=NULL;}
    if (RHS!=NULL)  { delete RHS; RHS=NULL;}
    if (DD!=NULL)   { delete DD; DD=NULL;}
    if (ND!=NULL)   { delete ND; ND=NULL;}
  }
  
  virtual std::string GetName() const=0;

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
  
  virtual void BasicInit(const std::string& paramfile) {
    _paramfile = paramfile;

    ConstructEquation         ();
    ConstructExactSolution    ();
    ConstructRightHandSideData();
    ConstructDirichletData    ();
    ConstructNeumannData      ();
    ConstructInitialCondition ();
    ConstructBoundaryManager();
  }

  const RightHandSideData*  GetRightHandSideData() const { return  RHS;}
  const DirichletData*      GetDirichletData()     const { return  DD;}
  const NeumannData*        GetNeumannData()       const { return  ND;}
  const InitialCondition*   GetInitialCondition() const { return IC;}
  const ExactSolution*   GetExactSolution   () const { return ES;}
  const Equation*        GetEquation        () const { return EQ;}
  const BoundaryManager* GetBoundaryManager () const { return BM;}

  void SetTime(double time, double dt) const {
    if (EQ) EQ -> SetTime(dt);
    if (ES) ES -> SetTime(time);
  }
};


#endif
