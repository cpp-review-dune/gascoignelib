#include  "stokes.h"
#include  "filescanner.h"

/*-----------------------------------------*/

Stokes::~Stokes()
{
}

/*-----------------------------------------*/

Stokes::Stokes() : Equation()
{
  penalty = 0.; visc = 1.;
}
 
/*-----------------------------------------*/

void Stokes::SetTimePattern(TimePattern& P) const
{
  P.reservesize(ncomp(),ncomp(),0.);
  P(0,0) = penalty;
  P(1,1) = 1.;
  P(2,2) = 1.;
}

/*-----------------------------------------*/

Stokes::Stokes(const std::string& filename) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &visc , 0.01);
  DFH.insert("penalty",&penalty, 0.);

  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

double Stokes::Laplace(const DerivativeVector& U, 
			       const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y();
}

/*-----------------------------------------*/

double Stokes::Divergence(const FemFunction& U) const
{
  return U[1].x() + U[2].y();
}

/*-----------------------------------------*/

void Stokes::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  b[1] -= U[0].m()*N.x();
  b[2] -= U[0].m()*N.y();
	  
  // viscous terms
  b[1] += visc * Laplace(U[1],N);
  b[2] += visc * Laplace(U[2],N);
}

/*-----------------------------------------*/

void Stokes::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  A(0,1) += M.x()*N.m();
  A(0,2) += M.y()*N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();

  double laplace = Laplace(M,N);
  A(1,1) += visc*laplace;
  A(2,2) += visc*laplace;
}

/*------------------------------------------------------*/
