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
public:


private:


protected:

  mutable double _us, _vs;

public:


//
///  Constructor 
//
  LocalInitialCondition(const LocalEquation* EQ);
  
  std::string GetName() const { return "Local";}
  int GetNcomp() const {return 2;}  
  double operator()(int c, const Vertex2d& v) const;
};


#endif
