#ifndef  __ProblemDescriptorBackward_h
#define  __ProblemDescriptorBackward_h

#include  "problemdescriptorinterface.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////




class ProblemDescriptorBackward : public ProblemDescriptorInterface
{
public:


private:


protected:

  void ConstructEquation();
  
public:


//
///  Constructor 
//
    ProblemDescriptorBackward() {}
    
    std::string GetName() const {return "Backward";}

};


#endif
