#ifndef  __ZeroDirichletData_h
#define  __ZeroDirichletData_h

#include  "dirichletdata.h"

/*-----------------------------------------*/


class ZeroDirichletData : public DirichletData
{
protected:


public:

  ZeroDirichletData() : DirichletData() {}
  std::string GetName() const {return "Zero";}
  void operator()(Gascoigne::DoubleVector& b, const Vertex2d& v, int col) const {b.zero();}
};


#endif
