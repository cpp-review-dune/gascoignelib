#include  "laplace.h"
#include  "filescanner.h"

using namespace Gascoigne;

/*-----------------------------------------*/

Laplace:: Laplace()
{
  gamma = 0.;
  visc = 1.;
}

/*-----------------------------------------*/

Laplace::Laplace(const ParamFile* pf) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("gamma",&gamma,0.);
  DFH.insert("visc",&visc,1.);
  FileScanner FS(DFH,pf,"Equation");
}
 
/*-----------------------------------------*/

void Laplace::OperatorStrong(Vector& b, const FemFunction& U)const
{
  b[0] -= visc*U[0].D();
}
 
/*-----------------------------------------*/

void Laplace::SetTimePattern(TimePattern& P) const
{
  P.reservesize(ncomp(),ncomp(),0.);
  P(0,0) = 1.;
}

/*-----------------------------------------*/

void Laplace::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  b[0] += visc* (U[0].x()*N.x()+U[0].y()*N.y());
}

/*-----------------------------------------*/

void Laplace::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  A(0,0) += visc* (M.x()*N.x()+M.y()*N.y());
}

/*-----------------------------------------*/
