#include  "curvedboundary.h"

/*---------------------------------------------*/

double SinusBoundary::operator()(const Vertex2d& c) const 
{
  double pi = 3.1415925359;
  return -1.+0.25*sin(pi/2.*(c.x()+1.)) - c.y();
}
