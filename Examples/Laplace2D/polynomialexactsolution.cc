#include  "polynomialexactsolution.h"


using namespace std;

/*-----------------------------------------*/

double PolynomialExactSolution::operator()(int c, const Vertex2d& v)const 
{
//   return v.x()*(1.-v.x())*v.y()*(1.-v.y());
  return v.x()*(1.-v.x())*v.y()*(1.-v.y()) *  (exp(v.x()+v.y()));
}
