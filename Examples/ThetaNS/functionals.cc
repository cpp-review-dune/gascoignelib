#include "functionals.h"


#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"

extern double __DT;
extern double __TIME;

namespace Gascoigne
{




  void LiftBRHS::operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v, const Vertex2d& n, int color) const 
  {
    if (__TIME<=5.0) return;
    if (__TIME>=6.0) return;
    if (color!=80) return;
    b[0] += 20.0 * N.m() * n.y();
    b[2] -= 20.0 * __nu * (n.x()*N.x() + n.y() * N.y());
  }




}
