#ifndef  __boundaryfunction_h
#define  __boundaryfunction_h

#include  "numfixarray.h"

/*---------------------------------------------------*/

template<int DIM>
class BoundaryFunction
{
protected:
  
  typedef numfixarray<DIM,double>  Vector;
  
public :
  
  virtual ~BoundaryFunction() {}
  virtual double operator()(const Vector& c) const=0;
  
  virtual void grad(Vector& dst, const Vector& src) const;
  void newton(Vector&) const;
};

/*---------------------------------------------------*/

#endif
