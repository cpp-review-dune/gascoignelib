#ifndef  __LocalInitialCondition_h
#define  __LocalInitialCondition_h



/////////////////////////////////////////////
///
///@brief
///  ... comments LocalInitialCondition

///
///
/////////////////////////////////////////////


#include  "domainrighthandside.h"
#include  "localequation.h"

class LocalInitialCondition : public Gascoigne::DomainRightHandSide
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
  double operator()(int c, const Gascoigne::Vertex2d& v) const;
};


#endif
