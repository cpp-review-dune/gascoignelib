#ifndef  __ProblemDescriptorInitial_h
#define  __ProblemDescriptorInitial_h

#include  "problemdescriptorinterface.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////




class ProblemDescriptorInitial : public ProblemDescriptorInterface
{
public:


private:


protected:

  void ConstructEquation();
  void ConstructInitialCondition();
  
public:


//
///  Constructor 
//
    ProblemDescriptorInitial() {}
    
    std::string GetName() const {return "Initial";}

};


#endif
