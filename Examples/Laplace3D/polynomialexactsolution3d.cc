#include  "polynomialexactsolution3d.h"

using namespace std;

/*-----------------------------------------*/

double PolynomialExactSolution3d::operator()(int c, const Vertex3d& v)const 
{
  return quadratic(v.x()) * quadratic(v.y()) * quadratic(v.z());
}
