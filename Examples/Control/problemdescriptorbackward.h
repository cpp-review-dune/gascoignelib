#ifndef  __ProblemDescriptorBackward_h
#define  __ProblemDescriptorBackward_h

#include  "problemdescriptorbase.h"
#include  "backwardequation.h"

/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////


class ProblemDescriptorBackward : public Gascoigne::ProblemDescriptorBase
{
 public:

  std::string GetName() const {return "Backward";}

  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetEquationPointer() = new BackwardEquation(GetParamFile());
    ProblemDescriptorBase::BasicInit(pf);
  }

};


#endif
