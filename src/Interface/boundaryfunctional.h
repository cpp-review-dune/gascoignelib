#ifndef  __BoundaryFunctional_h
#define  __BoundaryFunctional_h

#include  "nvector.h"
#include  <string>
#include  <set>
#include  "functional.h"

/*-----------------------------------------*/


class BoundaryFunctional : public virtual Functional
{

public:

  BoundaryFunctional() : Functional() {}

  virtual ~BoundaryFunctional() {};

  virtual std::set<int> GetColors() const=0;
  virtual double J(const Gascoigne::FemFunction& U, const Vertex2d& v) const {assert(0);}
  virtual void J(Vector& b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const {assert(0);}

};


#endif
