#include  "stokes3d.h"

using namespace Gascoigne;

/*-----------------------------------------*/

Stokes3d::~Stokes3d()
{
}

/*-----------------------------------------*/

Stokes3d::Stokes3d() : Stokes2d()
{}

/*-----------------------------------------*/

Stokes3d::Stokes3d(const ParamFile* pf) : Stokes2d(pf)
{
}

/*-----------------------------------------*/

double Stokes3d::Divergence(const FemFunction& U) const
{
  return U[1].x() + U[2].y() + U[3].z();
}

/*-----------------------------------------*/

double Stokes3d::Laplace(const TestFunction& U, const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y() + U.z()*N.z();
}

/*-----------------------------------------*/

void Stokes3d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity //////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////

  b[1] -= U[0].m()*N.x();
  b[2] -= U[0].m()*N.y();
  b[3] -= U[0].m()*N.z();
	  
  // viscous terms
  b[1] += visc * Laplace(U[1],N);
  b[2] += visc * Laplace(U[2],N);
  b[3] += visc * Laplace(U[3],N);
}

/*-----------------------------------------*/

void Stokes3d::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  double laplace = Laplace(M,N);
     
  ////////////// Continuity /////////////////////////////

  A(0,0) += M.m()*N.m() * penalty;
  A(0,1) += M.x()*N.m();
  A(0,2) += M.y()*N.m();
  A(0,3) += M.z()*N.m();

  ////////////// Momentum ///////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();
  A(3,0) -= M.m()*N.z();

  A(1,1) += visc*laplace;
  A(2,2) += visc*laplace;
  A(3,3) += visc*laplace;
}

/*-----------------------------------------*/
