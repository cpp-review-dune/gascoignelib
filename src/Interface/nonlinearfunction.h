#ifndef  __NonlinearFunction_h
#define  __NonlinearFunction_h

#include  "nvector.h"
#include  "vertex.h"

/*-----------------------------------------*/


class NonlinearFunction
{
protected:

  typedef  nvector<double>           Vector;

public:

  virtual ~NonlinearFunction() {};

  virtual double operator()(const Vertex2d& v, const FemFunction& U) const=0;
  virtual double operator()(const Vertex3d& v, const FemFunction& U) const {};
};

#endif
