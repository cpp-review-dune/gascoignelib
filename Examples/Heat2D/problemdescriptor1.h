#ifndef  __ProblemDescriptor1_h
#define  __ProblemDescriptor1_h

#include  "problemdescriptorinterface.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor1

///
///
/////////////////////////////////////////////




class ProblemDescriptor1 : public ProblemDescriptorInterface
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
    ProblemDescriptor1() {}
    
    std::string GetName() const {return "Local";}

};


#endif
