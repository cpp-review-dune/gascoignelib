#include  "stokesgls3d.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
StokesGls3d::~StokesGls3d()
{
}

/*-----------------------------------------*/

StokesGls3d::StokesGls3d() : Stokes3d(), GlsEquation()
{
}
 
/*-----------------------------------------*/

StokesGls3d::StokesGls3d(const ParamFile* pf) : Stokes3d(), GlsEquation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH, pf, "Equation");
}

/*-----------------------------------------*/

void StokesGls3d::glspoint(double h, const FemFunction& U, const Vertex3d& v)const
{
  ST.ReInit(h,_visc);
}

/*-----------------------------------------*/

void StokesGls3d::L(DoubleVector& dst, const FemFunction& U) const
{
  dst[0] = Divergence(U);
  dst[1] = U[0].x();
  dst[2] = U[0].y();
  dst[3] = U[0].z();
}

/*-----------------------------------------*/

void StokesGls3d::S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const
{
  dst(0,1) = ST.alpha() * N.x();
  dst(0,2) = ST.alpha() * N.y();
  dst(0,3) = ST.alpha() * N.z();
}

/*-----------------------------------------*/

void StokesGls3d::LMatrix(nmatrix<double>& A, const FemFunction& U, const TestFunction& N) const
{
  A(0,1) = N.x();
  A(0,2) = N.y();
  A(0,3) = N.z();
  A(1,0) = N.x();
  A(2,0) = N.y();
  A(3,0) = N.z();
}
}

/*-----------------------------------------*/

