#ifndef  __ProblemDescriptorForward_h
#define  __ProblemDescriptorForward_h

#include  "problemdescriptorinterface.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////




class ProblemDescriptorForward : public ProblemDescriptorInterface
{
public:


private:


protected:

  void ConstructEquation();
  
public:


//
///  Constructor 
//
    ProblemDescriptorForward() {}
    
    std::string GetName() const {return "Forward";}

};


#endif
