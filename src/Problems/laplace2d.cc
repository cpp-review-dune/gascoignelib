#include  "laplace2d.h"
#include  "filescanner.h"

using namespace Gascoigne;

/*-----------------------------------------*/

Laplace2d:: Laplace2d()
{
  gamma = 0.;
  visc = 1.;
}

/*-----------------------------------------*/

Laplace2d::Laplace2d(const ParamFile* pf) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("gamma",&gamma,0.);
  DFH.insert("visc",&visc,1.);
  FileScanner FS(DFH,pf,"Equation");
}
 
/*-----------------------------------------*/

void Laplace2d::OperatorStrong(DoubleVector& b, const FemFunction& U)const
{
  b[0] -= visc*U[0].D();
}
 
/*-----------------------------------------*/

void Laplace2d::SetTimePattern(TimePattern& P) const
{
  P.reservesize(ncomp(),ncomp(),0.);
  P(0,0) = 1.;
}

/*-----------------------------------------*/

void Laplace2d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  b[0] += visc* (U[0].x()*N.x()+U[0].y()*N.y());
}

/*-----------------------------------------*/

void Laplace2d::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  A(0,0) += visc* (M.x()*N.x()+M.y()*N.y());
}

/*-----------------------------------------*/
