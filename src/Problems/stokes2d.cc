#include  "stokes2d.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
Stokes2d::~Stokes2d()
{
}

/*-----------------------------------------*/

Stokes2d::Stokes2d() : Equation()
{
  _penalty = 0.; _visc = 1.;
}
 
/*-----------------------------------------*/

void Stokes2d::SetTimePattern(TimePattern& P) const
{
  P.reservesize(GetNcomp(),GetNcomp(),0.);
  P(0,0) = _penalty;
  P(1,1) = 1.;
  P(2,2) = 1.;
}

/*-----------------------------------------*/

Stokes2d::Stokes2d(const ParamFile* pf) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 1.);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

double Stokes2d::Laplace(const TestFunction& U, 
			       const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y();
}

/*-----------------------------------------*/

double Stokes2d::Divergence(const FemFunction& U) const
{
  return U[1].x() + U[2].y();
}

/*-----------------------------------------*/

void Stokes2d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  b[1] -= U[0].m()*N.x();
  b[2] -= U[0].m()*N.y();
	  
  // viscous terms
  b[1] += _visc * Laplace(U[1],N);
  b[2] += _visc * Laplace(U[2],N);
}

/*-----------------------------------------*/

void Stokes2d::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  A(0,1) += M.x()*N.m();
  A(0,2) += M.y()*N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();

  double laplace = Laplace(M,N);
  A(1,1) += _visc*laplace;
  A(2,2) += _visc*laplace;
}
}

/*------------------------------------------------------*/
