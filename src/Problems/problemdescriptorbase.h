#ifndef __ProblemDescriptorBase_h
#define __ProblemDescriptorBase_h

#include "problemdescriptorinterface.h"

/*--------------------------------------------------*/

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
  
  const ParamFile*& GetParamFilePointer()           { return _paramfile;}
  
  Equation*&          GetEquationPointer()          { return EQ;}
  BoundaryManager*&   GetBoundaryManagerPointer()   { return BM;}
  ExactSolution*&     GetExactSolutionPointer()     { return ES;}
  Application*&       GetInitialConditionPointer()  { return IC;}
  Application*&       GetRightHandSidePointer()     { return RHS;}
  DirichletData*&     GetDirichletDataPointer()     { return DD;}
  BoundaryRightHandSide*& GetBoundaryRightHandSidePointer() { return BRHS;}
  BoundaryEquation*&      GetBoundaryEquationPointer()      { return BE; }

  BoundaryManager* GetBoundaryManager () { return BM;}
 public:

  ProblemDescriptorBase();

  ~ProblemDescriptorBase();

  std::ostream& OutputSettings(std::ostream& os) const;
  
  void BasicInit(const ParamFile* pf);
  
  const ParamFile* GetParamFile() const {return _paramfile;}

  const Application*           GetRightHandSide        () const { return  RHS;}
  const DirichletData*         GetDirichletData        () const { return  DD;}
  const BoundaryRightHandSide* GetBoundaryRightHandSide() const { return  BRHS;}
  const BoundaryEquation*      GetBoundaryEquation     () const { return  BE;}
  const Application*           GetInitialCondition     () const { return IC;}
  const ExactSolution*         GetExactSolution        () const { return ES;}
  const Equation*              GetEquation             () const { return EQ;}
  const BoundaryManager*       GetBoundaryManager      () const { return BM;}

  void SetTime(double time, double dt) const;
};
}

/*--------------------------------------------------*/

#endif
