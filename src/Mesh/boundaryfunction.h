#ifndef  __boundaryfunction_h
#define  __boundaryfunction_h

#include  <string>
#include  "vertex.h"

/*---------------------------------------------------*/

template<int DIM>
class BoundaryFunction
{
protected:
  
  typedef Vertex<DIM>  Vector;
  
public :
  
  BoundaryFunction() {}
  virtual ~BoundaryFunction() {}

  virtual std::string GetName() const=0;
  virtual double operator()(const Vector& c) const=0;
  
  virtual void grad(Vector& dst, const Vector& src) const;
  void newton(Vector&) const;
};

/*---------------------------------------------------*/

#endif
