#include  "stokesgls.h"
#include  "filescanner.h"

using namespace Gascoigne;

/*-----------------------------------------*/

StokesGls::~StokesGls()
{
}

/*-----------------------------------------*/

StokesGls::StokesGls() : Stokes(), GlsEquation()
{}
 
/*-----------------------------------------*/

StokesGls::StokesGls(const ParamFile* pf) : Stokes(), GlsEquation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &visc , 1.);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("penalty",&penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

void StokesGls::glspoint
(double h, const FemFunction& U, const Vertex2d& v)const
{
  ST.ReInit(h,visc);
}

/*-----------------------------------------*/

void StokesGls::L(nvector<double>& dst, const FemFunction& U) const
{
  dst[0] = Divergence(U);
  dst[1] = U[0].x();
  dst[2] = U[0].y();
}

/*-----------------------------------------*/

void StokesGls::S(nmatrix<double>& dst, const FemFunction& U, 
		  const TestFunction& N) const
{
  dst(0,1) = ST.alpha() * N.x();
  dst(0,2) = ST.alpha() * N.y();
  // div-div
  //  dst(1,0) = ST.alpha() * N.x();
  //  dst(2,0) = ST.alpha() * N.y();
}

/*-----------------------------------------*/

void StokesGls::LMatrix(nmatrix<double>& A, 
			const FemFunction& U,
			const DerivativeVector& V) const
{
  A(0,1) = V.x();
  A(0,2) = V.y();
  A(1,0) = V.x();
  A(2,0) = V.y();
}

/*-----------------------------------------*/

