#include  "navierstokes.h"
#include  "filescanner.h"

using namespace Gascoigne;

/*-----------------------------------------*/

NavierStokes::~NavierStokes()
{
}

/*-----------------------------------------*/

NavierStokes::NavierStokes() : Equation()
{
  penalty = 0.; visc = 0.01; _h = 0.;
}

/*-----------------------------------------*/

NavierStokes::NavierStokes(const ParamFile* pf) : Equation()
{
  _h = 0.;
  DataFormatHandler DFH;
  DFH.insert("visc" , &visc , 1.);
  DFH.insert("cut"  , &cut  ,  1.e10);
  DFH.insert("penalty",&penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

void NavierStokes::OperatorStrong(Vector& b, const FemFunction& U) const
{
  b[0] = Divergence(U);
  b[1] = Convection(U,U[1]) - visc * U[1].D() + U[0].x();
  b[2] = Convection(U,U[2]) - visc * U[2].D() + U[0].y();
}

/*-----------------------------------------*/

double NavierStokes::Laplace(const DerivativeVector& U, 
				     const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y();
}

/*-----------------------------------------*/

double NavierStokes::Convection(const FemFunction& U, 
					const TestFunction& N) const
{
  return U[1].m()*N.x() + U[2].m()*N.y();
}

/*-----------------------------------------*/

double NavierStokes::Divergence(const FemFunction& U) const
{
  return U[1].x() + U[2].y();
}
 
/*-----------------------------------------*/

void NavierStokes::SetTimePattern(TimePattern& P) const
{
  P.reservesize(ncomp(),ncomp(),0.);
  P(0,0) = penalty;
  P(1,1) = 1.;
  P(2,2) = 1.;
}

/*-----------------------------------------*/

void NavierStokes::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  b[1] -= U[0].m()*N.x();
  b[2] -= U[0].m()*N.y();
	  
  b[1] += Convection(U,U[1]) * N.m();
  b[2] += Convection(U,U[2]) * N.m();

  // viscous terms
  b[1] += visc * Laplace(U[1],N);
  b[2] += visc * Laplace(U[2],N);
}

/*-----------------------------------------*/

void NavierStokes::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  double MN = M.m()*N.m();
  double Mx = M.x()*N.m();
  double My = M.y()*N.m();
  double laplace = Laplace(M,N);
     
  ////////////// Continuity ////////////////////////////////////////////////

//   A(0,0) += MN * alpha00;// * ST->h()* ST->h();
  A(0,1) += Mx;
  A(0,2) += My;

  ////////////// Momentum ////////////////////////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();

  double cl = Convection(U,M) * N.m() + visc*laplace;

  A(1,1) += cl;
  A(2,2) += cl;

  double Cut = cut * _h;

  A(1,1) += GascoigneMath::max(U[1].x()*MN, -Cut);
  A(2,2) += GascoigneMath::max(U[2].y()*MN, -Cut);

  A(1,2) += GascoigneMath::min(U[1].y()*MN, Cut);
  A(2,1) += GascoigneMath::min(U[2].x()*MN, Cut);
}


/*-----------------------------------------*/


