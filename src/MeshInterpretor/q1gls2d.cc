#include  "q1gls2d.h"
#include  "galerkinglsintegrator.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq12d.h"

using namespace std;

/* ----------------------------------------- */

void Q1Gls2d::BasicInit(const std::string& paramfile)
{
  assert(CellMeshInterpretor::GetIntegrator()==NULL);
  BasicMeshInterpretor::GetIntegratorPointer() =  new GalerkinGlsIntegrator<2>;

  assert(CellMeshInterpretor::GetFem()==NULL);
  typedef Transformation2d<BaseQ12d>           TransQ1;
  typedef FiniteElement<2,1,TransQ1,BaseQ12d>  FiniteElement;
  BasicMeshInterpretor::GetFemPointer() =  new FiniteElement;

  CellMeshInterpretor::BasicInit(paramfile);
}

/* ----------------------------------------- */
