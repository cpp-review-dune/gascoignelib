#ifndef  __ProblemDescriptorTerminal_h
#define  __ProblemDescriptorTerminal_h

#include  "problemdescriptorinterface.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////




class ProblemDescriptorTerminal : public ProblemDescriptorInterface
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
    ProblemDescriptorTerminal() {}
    
    std::string GetName() const {return "Terminal";}

};


#endif
