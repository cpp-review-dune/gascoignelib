#ifndef  __ProblemDescriptor_h
#define  __ProblemDescriptor_h

#include  "problemdescriptorinterface.h"

/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////

class ProblemDescriptor : public ProblemDescriptorInterface
{
public:


private:


protected:

  void ConstructEquation();
  void ConstructDirichletData();
  
public:

//
///  Constructor 
//
    ProblemDescriptor() {}
    
    std::string GetName() const {return "Local";}
};

#endif
