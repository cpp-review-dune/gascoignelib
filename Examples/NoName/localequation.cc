#include  "localequation.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

void LocalEquation::point(double h, const FemFunction& U, const Vertex2d& v) const 
{
  NavierStokes::point(h,U,v);
}

/*-----------------------------------------*/

void LocalEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  NavierStokes::Form(b,U,N);
}

/*-----------------------------------------*/

void LocalEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  NavierStokes::Matrix(A,U,M,N);
}

