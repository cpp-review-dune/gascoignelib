#ifndef  __ProblemDescriptorInitial_h
#define  __ProblemDescriptorInitial_h

#include  "problemdescriptorbase.h"
#include  "localequation.h"
#include  "localinitialcondition.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////


class ProblemDescriptorInitial : public ProblemDescriptorBase
{
public:
    std::string GetName() const {return "Initial";}

  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetEquationPointer() = new LocalEquation(GetParamFile());
    const LocalEquation* LEQ = dynamic_cast<const LocalEquation*>(GetEquation());
    GetInitialConditionPointer() = new LocalInitialCondition(LEQ);
    ProblemDescriptorBase::BasicInit(pf);
  }

};


#endif
