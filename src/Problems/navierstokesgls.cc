#include  "navierstokesgls.h"
#include  "filescanner.h"

using namespace std;

/*-----------------------------------------*/

NavierStokesGls::~NavierStokesGls()
{
}

/*-----------------------------------------*/

NavierStokesGls::NavierStokesGls() 
  : NavierStokes(), GlsEquation()
{
  penalty = 0.; visc = 0.01;
  
  ST.xeta0 = 6.;
  ST.delta0 = ST.alpha0 = 0.2;
  ST.xeta0 = 6.; 
}

/*-----------------------------------------*/

NavierStokesGls::NavierStokesGls(const std::string& paramfile) 
  : NavierStokes(), GlsEquation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &visc , 0.01);
  DFH.insert("cut"  , &cut  ,  1.e8);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("delta", &ST.delta0, 0.25);
  DFH.insert("tau"  , &ST.tau0  , 0.);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&penalty, 0.);

  FileScanner FS(DFH,paramfile,"Equation");
}

/*-----------------------------------------*/

void NavierStokesGls::glspoint(double h, const FemFunction& U, const Vertex2d& v)const
{
  ST.ReInit(h,visc,U[1].m(),U[2].m());
}

/*-----------------------------------------*/

void NavierStokesGls::L(nvector<double>& dst, const FemFunction& U) const
{
  dst[0] = Divergence(U);
  dst[1] = Convection(U,U[1]) + U[0].x();
  dst[2] = Convection(U,U[2]) + U[0].y();
}

/*-----------------------------------------*/

void NavierStokesGls::S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const
{
  dst(0,1) = ST.alpha() * N.x();
  dst(0,2) = ST.alpha() * N.y();
  dst(1,0) = ST.tau()   * N.x();
  dst(2,0) = ST.tau()   * N.y();

  dst(1,1) = ST.delta() * Convection(U,N);
  dst(2,2) = ST.delta() * Convection(U,N);
}

/*-----------------------------------------*/

void NavierStokesGls::LMatrix(nmatrix<double>& A, const FemFunction& U, const TestFunction& N) const
{
  A(0,1) = N.x();
  A(0,2) = N.y();
     
  A(1,0) = N.x();
  A(2,0) = N.y();

  double cl = Convection(U,N);

  A(1,1) = cl + U[1].x()*N.m();
  A(2,2) = cl + U[2].y()*N.m();
  A(1,2) =      U[1].y()*N.m();
  A(2,1) =      U[2].x()*N.m();
}

/*-----------------------------------------*/

void NavierStokesGls::SMatrix(nvector<double>& dst, const FemFunction& U, const FemFunction& M, const FemFunction& N) const
{
  dst[1] = ST.delta() * Convection(M,N[1]);
  dst[2] = ST.delta() * Convection(M,N[2]);
}
