#include  "neumanndatabyexactsolution.h"

using namespace std;

/* ----------------------------------------- */

void NeumannDataByExactSolution::operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v, const Vertex2d& n, int col) const
{
  b[0] += ( _ES->x(0,v)*n.x()+_ES->y(0,v)*n.y() ) * N.m();
}
