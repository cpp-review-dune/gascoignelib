#include "stokeslps3d.h"
#include "filescanner.h"

/*-----------------------------------------*/

namespace Gascoigne
{
StokesLps3d::~StokesLps3d() {}

/*-----------------------------------------*/

StokesLps3d::StokesLps3d() : 
  Stokes3d(), LpsEquation()
{ 
  _penalty = 0.; _visc = 1.;
  
  ST.xeta0 = 6.;
  ST.alpha0 = 0.2;
}

/*-----------------------------------------*/

StokesLps3d::StokesLps3d(const ParamFile* filename) : 
  Stokes3d(), LpsEquation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

void StokesLps3d::lpspoint(double h, const FemFunction& U, const Vertex3d& v)const
{
  ST.ReInit(h,_visc);
}

/*-----------------------------------------*/

void StokesLps3d::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
{
  b[0] += ST.alpha() * Laplace(UP[0],N);
}
 
/*-----------------------------------------*/

void StokesLps3d::StabMatrix(EntryMatrix& A,  const FemFunction& U, const TestFunction& Np, 
			     const TestFunction& Mp) const
{
  A(0,0) += ST.alpha() * Laplace(Mp,Np);
}
}
