#ifndef  __OneDirichletData_h
#define  __OneDirichletData_h

#include  "dirichletdata.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class OneDirichletData : public DirichletData
{
protected:


public:

  OneDirichletData() : DirichletData() {}
  std::string GetName() const {return "One";}
  void operator()(Vector& b, const Vertex2d& v, int col) const {b.zero(); b[0] = 1.;}
};
}

#endif
