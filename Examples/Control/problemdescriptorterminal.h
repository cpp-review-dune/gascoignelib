#ifndef  __ProblemDescriptorTerminal_h
#define  __ProblemDescriptorTerminal_h

#include  "problemdescriptorbase.h"
#include  "localequation.h"
#include  "localterminalcondition.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////




class ProblemDescriptorTerminal : public ProblemDescriptorBase
{
public:
    std::string GetName() const {return "Terminal";}

  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetEquationPointer() = new LocalEquation(GetParamFile());
    GetInitialConditionPointer() = new LocalTerminalCondition();
    ProblemDescriptorBase::BasicInit(pf);
  }

};


#endif
