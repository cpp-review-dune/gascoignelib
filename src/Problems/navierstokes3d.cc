#include  "navierstokes3d.h"

/*-----------------------------------------*/

NavierStokes3d::~NavierStokes3d() {}

/*-----------------------------------------*/

NavierStokes3d::NavierStokes3d() : NavierStokes() {}

/*-----------------------------------------*/

void NavierStokes3d::SetTimePattern(TimePattern& P) const
{
  P.reservesize(ncomp(),ncomp(),0.);
  P(0,0) = penalty;
  P(1,1) = 1.;
  P(2,2) = 1.;
  P(3,3) = 1.;
}

/*-----------------------------------------*/

NavierStokes3d::NavierStokes3d(const std::string& filename) 
  : NavierStokes(filename) {}


/*-----------------------------------------*/

void NavierStokes3d::OperatorStrong(Vector& b, const FemFunction& U) const
{
  b[0] = Divergence(U);
  b[1] = Convection(U,U[1]) - visc * U[1].D() + U[0].x();
  b[2] = Convection(U,U[2]) - visc * U[2].D() + U[0].y();
  b[3] = Convection(U,U[3]) - visc * U[3].D() + U[0].z();
}

/*-----------------------------------------*/

void NavierStokes3d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  b[1] -= U[0].m()*N.x();
  b[2] -= U[0].m()*N.y();
  b[3] -= U[0].m()*N.z();
	  
  b[1] += Convection(U,U[1]) * N.m();
  b[2] += Convection(U,U[2]) * N.m();
  b[3] += Convection(U,U[3]) * N.m();

  // viscous terms
  b[1] += visc * Laplace(U[1],N);
  b[2] += visc * Laplace(U[2],N);
  b[3] += visc * Laplace(U[3],N);
}

/*-----------------------------------------*/

void NavierStokes3d::Matrix
(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  double MN = M.m()*N.m();
     
  ////////////// Continuity ////////////////////////////////////////////////

  A(0,0) += MN * penalty;
  A(0,1) += M.x()*N.m();
  A(0,2) += M.y()*N.m();
  A(0,3) += M.z()*N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();
  A(3,0) -= M.m()*N.z();

  double sum = Convection(U,M)*N.m() + visc*Laplace(M,N);
  A(1,1) += sum;
  A(2,2) += sum;
  A(3,3) += sum;

  double tau = cut * _h;

  A(1,1) += GascoigneMath::max(U[1].x()*MN, -tau);
  A(2,2) += GascoigneMath::max(U[2].y()*MN, -tau);
  A(3,3) += GascoigneMath::max(U[3].z()*MN, -tau);

//   A(1,2) += GascoigneMath::min(U[1].y()*MN, tau);
//   A(1,3) += GascoigneMath::min(U[1].z()*MN, tau);
//   A(2,1) += GascoigneMath::min(U[2].x()*MN, tau);
//   A(2,3) += GascoigneMath::min(U[2].z()*MN, tau);
//   A(3,1) += GascoigneMath::min(U[3].x()*MN, tau);
//   A(3,2) += GascoigneMath::min(U[3].y()*MN, tau);
}

/*-----------------------------------------*/

std::string NavierStokes3d::GetName() const 
{ 
  return "NavierStokes3d";
}

/*-----------------------------------------*/
