#include  "stokesgls2d.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
StokesGls2d::~StokesGls2d()
{
}

/*-----------------------------------------*/

StokesGls2d::StokesGls2d() : Stokes2d(), GlsEquation()
{
  _penalty = 0.; 
  _visc = 1.;
  ST.alpha0 = 1.;
}
 
/*-----------------------------------------*/

StokesGls2d::StokesGls2d(const ParamFile* pf) : Stokes2d(), GlsEquation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 1.);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

void StokesGls2d::glspoint(double h, const FemFunction& U, const Vertex2d& v)const
{
  ST.ReInit(h,_visc);
}

/*-----------------------------------------*/

void StokesGls2d::L(DoubleVector& dst, const FemFunction& U) const
{
  dst[0] = Divergence(U);
  dst[1] = U[0].x();
  dst[2] = U[0].y();
}

/*-----------------------------------------*/

void StokesGls2d::S(nmatrix<double>& dst, const FemFunction& U, 
		  const TestFunction& N) const
{
  dst(0,1) = ST.alpha() * N.x();
  dst(0,2) = ST.alpha() * N.y();
  // div-div
  //  dst(1,0) = ST.alpha() * N.x();
  //  dst(2,0) = ST.alpha() * N.y();
}

/*-----------------------------------------*/

void StokesGls2d::LMatrix(nmatrix<double>& A, 
			const FemFunction& U,
			const TestFunction& V) const
{
  A(0,1) = V.x();
  A(0,2) = V.y();
  A(1,0) = V.x();
  A(2,0) = V.y();
}
}

/*-----------------------------------------*/

