#include  "navierstokesgls3d.h"
#include  "filescanner.h"

using namespace Gascoigne;

/*-----------------------------------------*/

NavierStokesGls3d::~NavierStokesGls3d() {}

/*-----------------------------------------*/

NavierStokesGls3d::NavierStokesGls3d() : NavierStokes3d(), GlsEquation() {}

/*-----------------------------------------*/

NavierStokesGls3d::NavierStokesGls3d(const ParamFile* pf) 
  : NavierStokes3d(), GlsEquation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &visc , 0.01);
  DFH.insert("cut"  , &cut  ,  1.e8);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("delta", &ST.delta0, 0.25);
  DFH.insert("tau"  , &ST.tau0  , 0.);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

std::string NavierStokesGls3d::GetName() const
{
  return "NavierStokesGls3d";
}

/*-----------------------------------------*/

void NavierStokesGls3d::glspointmatrix
(double h, const FemFunction& U, FemData& Q, const Vertex3d& v) const
{
  _h = h;
  ST.ReInit(h,visc,U[1].m(),U[2].m(),U[3].m());
}

/*-----------------------------------------*/

void NavierStokesGls3d::glspoint
(double h, const FemFunction& U, const Vertex3d& v) const
{
  _h = h;
  ST.ReInit(h,visc,U[1].m(),U[2].m(),U[3].m());
}
  
/*-----------------------------------------*/

void NavierStokesGls3d::L(nvector<double>& dst, const FemFunction& U) const
{
  dst[0] = Divergence(U);
  dst[1] = Convection(U,U[1]) + U[0].x();
  dst[2] = Convection(U,U[2]) + U[0].y();
  dst[3] = Convection(U,U[3]) + U[0].z();
}

/*-----------------------------------------*/

void NavierStokesGls3d::S(nmatrix<double>& dst, const FemFunction& U, 
			     const TestFunction& N) const
{
  dst(0,1) = ST.alpha() * N.x();
  dst(0,2) = ST.alpha() * N.y();
  dst(0,3) = ST.alpha() * N.z();
  dst(1,0) = ST.tau() * N.x();
  dst(2,0) = ST.tau() * N.y();
  dst(3,0) = ST.tau() * N.z();

  double conv = ST.delta() * Convection(U,N);
  
  dst(1,1) = conv;
  dst(2,2) = conv;
  dst(3,3) = conv;
}

/*-----------------------------------------*/

void NavierStokesGls3d::LMatrix(nmatrix<double>& A, 
				   const FemFunction& U,
				   const TestFunction& N) const
{
  A(0,1) = N.x();
  A(0,2) = N.y();
  A(0,3) = N.z();
     
  A(1,0) = N.x();
  A(2,0) = N.y();
  A(3,0) = N.z();

  double cl = Convection(U,N);

  A(1,1) = cl + U[1].x()*N.m();
  A(2,2) = cl + U[2].y()*N.m();
  A(3,3) = cl + U[3].z()*N.m();
  A(1,2) =      U[1].y()*N.m();
  A(1,3) =      U[1].z()*N.m();
  A(2,1) =      U[2].x()*N.m();
  A(2,3) =      U[2].z()*N.m();
  A(3,1) =      U[3].x()*N.m();
  A(3,2) =      U[3].y()*N.m();
}

/*-----------------------------------------*/
