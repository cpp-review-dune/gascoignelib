#ifndef  __FunctionalZero_h
#define  __FunctionalZero_h

#include  "functional.h"

/*-----------------------------------------*/


class ZeroFunctional : public Functional
{
protected:
 
public:


  ZeroFunctional() : Functional() {}
  ~ZeroFunctional() {}

  std::string GetName() const {return "zero";}

  double J(const Gascoigne::FemFunction& U, const Vertex2d& v) const{return 0.;}
};


#endif
