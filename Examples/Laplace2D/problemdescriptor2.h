#ifndef  __ProblemDescriptor2_h
#define  __ProblemDescriptor2_h

#include  "problemdescriptorinterface.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments ProblemDescriptor

///
///
/////////////////////////////////////////////




class ProblemDescriptor2 : public ProblemDescriptorInterface
{
public:


private:


protected:

  void ConstructEquation();
  void ConstructRightHandSideData();
  void ConstructDirichletData();
  void ConstructNeumannData();
  void ConstructExactSolution();
  
public:


//
///  Constructor 
//
    ProblemDescriptor2() {}
    
    std::string GetName() const {return "Local";}

};


#endif
