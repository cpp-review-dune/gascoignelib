#include  "navierstokeslps2d.h"
#include  "filescanner.h"

namespace Gascoigne
{

/*-----------------------------------------*/

NavierStokesLps2d::NavierStokesLps2d(const ParamFile* filename) 
  : LpsEquation(), NavierStokes2d() 
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("cut"  , &_cut  ,  1.e8);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("delta", &ST.delta0, 0.25);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

void NavierStokesLps2d::lpspoint(double h, const FemFunction& U, const Vertex2d& v)const
{
  ST.ReInit(h,_visc,U[1].m(),U[2].m());
}
 
/*-----------------------------------------*/

void NavierStokesLps2d::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
{
  b[0] += ST.alpha() * (UP[0].x()*N.x()+UP[0].y()*N.y());
  double betaN = Convection(U,N);
  double betaU1 = U[1].m()*UP[1].x()+U[2].m()*UP[1].y();
  double betaU2 = U[1].m()*UP[2].x()+U[2].m()*UP[2].y();
  b[1] += ST.delta() * betaU1*betaN;
  b[2] += ST.delta() * betaU2*betaN;
}
 
/*-----------------------------------------*/

void NavierStokesLps2d::StabMatrix(EntryMatrix& A,  const FemFunction& U, 
 const TestFunction& Np, const TestFunction& Mp) const
{
  double laplace = Laplace(Mp,Np);
  double betaM = Convection(U,Mp);
  double betaN = Convection(U,Np);
  double betabeta = ST.delta() * betaM*betaN;

  A(0,0) += ST.alpha() * laplace;
  A(1,1) += betabeta;
  A(2,2) += betabeta;
}
}
