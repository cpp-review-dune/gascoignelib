#ifndef  __DomainFunctional_h
#define  __DomainFunctional_h


/////////////////////////////////////////////
///
///@brief
///  ... comments DomainFunctional

///
///
/////////////////////////////////////////////

#include  "functional.h"

class DomainFunctional : public virtual Functional
{
public:

  DomainFunctional() : Functional() {}
  virtual ~DomainFunctional() {}

  virtual double J(const Gascoigne::FemFunction& U, const Vertex2d& v) const{ std::cerr << "DomainFunctional::J(Vertex2d&) not written\n"; abort(); return 0;}

  virtual double J(const Gascoigne::FemFunction& U, const Vertex3d& v) const{ std::cerr << "DomainFunctional::J(Vertex3d&) not written\n"; abort(); return 0;}

};


#endif
