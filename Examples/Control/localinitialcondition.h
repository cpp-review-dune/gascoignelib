#ifndef  __LocalInitialCondition_h
#define  __LocalInitialCondition_h



/////////////////////////////////////////////
///
///@brief
///  ... comments LocalInitialCondition

///
///
/////////////////////////////////////////////


#include  "initialcondition.h"
#include  "localequation.h"

class LocalInitialCondition : public InitialCondition
{
protected:

  int ncomp;

public:


//
///  Constructor 
//
  LocalInitialCondition(const LocalEquation* EQ);
  
  std::string GetName() const { return "Local";}
  int GetNcomp() const {return ncomp;}  
  double operator()(int c, const Vertex2d& v) const;
};


#endif