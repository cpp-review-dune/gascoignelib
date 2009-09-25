#ifndef __ProblemDescriptorBase_h
#define __ProblemDescriptorBase_h

#include "problemdescriptorinterface.h"

/*--------------------------------------------------*/

namespace Gascoigne
{
class ProblemDescriptorBase : public ProblemDescriptorInterface
{
 private:
  
  Equation                 *EQ;
  FaceEquation             *FEQ;
  BoundaryManager          *BM;
  ExactSolution            *ES;
  Application              *RHS;
  Application              *IC;
  DirichletData            *DD;
  PeriodicData             *PD;
  BoundaryRightHandSide    *BRHS;
  BoundaryInitialCondition *BIC;
  BoundaryEquation         *BE;
  ComponentInformation     *CI;

  
  const ParamFile *_paramfile;
  
 protected:
  
  const ParamFile*& GetParamFilePointer()           { return _paramfile;}
  
  Equation*&                 GetEquationPointer()                 { return EQ;}
  FaceEquation*&             GetFaceEquationPointer()             { return FEQ;}
  BoundaryManager*&          GetBoundaryManagerPointer()          { return BM;}
  ExactSolution*&            GetExactSolutionPointer()            { return ES;}
  Application*&              GetInitialConditionPointer()         { return IC;}
  Application*&              GetRightHandSidePointer()            { return RHS;}
  DirichletData*&            GetDirichletDataPointer()            { return DD;}
  PeriodicData*&             GetPeriodicDataPointer()             { return PD;}
  BoundaryRightHandSide*&    GetBoundaryRightHandSidePointer()    { return BRHS;}
  BoundaryInitialCondition*& GetBoundaryInitialConditionPointer() { return BIC;}
  BoundaryEquation*&         GetBoundaryEquationPointer()         { return BE; }

  BoundaryManager* GetBoundaryManager () { return BM;}

  ComponentInformation*&  GetComponentInformationPointer ()        { return CI;}

 public:

  ProblemDescriptorBase();

  ~ProblemDescriptorBase();

  std::ostream& OutputSettings(std::ostream& os) const;
  
  void BasicInit(const ParamFile* pf);
  
  const ParamFile* GetParamFile() const {return _paramfile;}

  const Application*              GetRightHandSide           () const { return RHS;}
  const DirichletData*            GetDirichletData           () const { return DD;}
  const PeriodicData*             GetPeriodicData            () const { return PD;}
  const BoundaryRightHandSide*    GetBoundaryRightHandSide   () const { return BRHS;}
  const BoundaryInitialCondition* GetBoundaryInitialCondition() const { return BIC;}
  const BoundaryEquation*         GetBoundaryEquation        () const { return BE;}
  const Application*              GetInitialCondition        () const { return IC;}
  const ExactSolution*            GetExactSolution           () const { return ES;}
  const Equation*                 GetEquation                () const { return EQ;}
  const FaceEquation*             GetFaceEquation            () const { return FEQ;}
  const BoundaryManager*          GetBoundaryManager         () const { return BM;}
  const ComponentInformation*     GetComponentInformation    () const { return CI;}

  void SetTime(double time, double dt) const;
};


 
}

/*--------------------------------------------------*/

#endif
