#ifndef  __ZeroDirichletData_h
#define  __ZeroDirichletData_h

#include  "dirichletdata.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class ZeroDirichletData : public DirichletData
{
protected:


public:

  ZeroDirichletData() : DirichletData() {}
  std::string GetName() const {return "Zero";}
  void operator()(DoubleVector& b, const Vertex2d& v, int col) const {b.zero();}
  void operator()(DoubleVector& b, const Vertex3d& v, int col) const {b.zero();}
};
}

#endif
