#include  "navierstokeslps3d.h"
#include  "filescanner.h"

namespace Gascoigne
{

/*-----------------------------------------*/

NavierStokesLps3d::NavierStokesLps3d(const ParamFile* filename) 
  : NavierStokes3d()
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

std::string NavierStokesLps3d::GetName() const
{
  return "NavierStokesLps3d";
}

/*-----------------------------------------*/

void NavierStokesLps3d::lpspoint
(double _h, const FemFunction& U, FemData& Q, const Vertex3d& v)const
{
  double h = _h;
  ST.ReInit(h,_visc,U[1].m(),U[2].m(),U[3].m());
}
 
/*-----------------------------------------*/

void NavierStokesLps3d::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
{
  b[0] += ST.alpha() * Laplace(UP[0],N);//(UP[0].x()*N.x()+UP[0].y()*N.y()+UP[0].z()*N.z());
  double betaN  = Convection(U,N);
  double betaU1 = Convection(U,UP[1]);//U[1].m()*UP[1].x()+U[2].m()*UP[1].y()+U[3].m()*UP[1].z();
  double betaU2 = Convection(U,UP[2]);//U[1].m()*UP[2].x()+U[2].m()*UP[2].y()+U[3].m()*UP[2].z();
  double betaU3 = Convection(U,UP[3]);//U[1].m()*UP[3].x()+U[2].m()*UP[3].y()+U[3].m()*UP[3].z();
  b[1] += ST.delta() * betaU1*betaN;
  b[2] += ST.delta() * betaU2*betaN;
  b[3] += ST.delta() * betaU3*betaN;
}
 
/*-----------------------------------------*/

void NavierStokesLps3d::StabMatrix(EntryMatrix& A,  const FemFunction& U, 
 const TestFunction& Np, const DerivativeVector& Mp) const
{
  double laplace = Laplace(Mp,Np);
  double betaM = Convection(U,Mp);
  double betaN = Convection(U,Np);
  double betabeta = ST.delta() * betaM*betaN;

  ////////////// Continuity ////////////////////////////////////////////////

  A(0,0) += ST.alpha() * laplace;
  A(1,1) += betabeta;
  A(2,2) += betabeta;
  A(3,3) += betabeta;
}
}

/*-----------------------------------------*/
