#ifndef  __ProblemDescriptorForward_h
#define  __ProblemDescriptorForward_h

#include  "problemdescriptorbase.h"
#include  "localequation.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////


class ProblemDescriptorForward : public Gascoigne::ProblemDescriptorBase
{
public:
    
  std::string GetName() const {return "Forward";}

  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new LocalEquation(GetParamFile());
    ProblemDescriptorBase::BasicInit(pf);
  }
  
};


#endif
