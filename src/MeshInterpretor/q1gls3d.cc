#include  "q1gls3d.h"
#include  "galerkinglsintegrator.h"
#include  "transformation3d.h"
#include  "finiteelement.h"
#include  "baseq13d.h"

using namespace std;

/* ----------------------------------------- */

void Q1Gls3d::BasicInit(const ParamFile* pf)
{
  assert(HN==NULL);
  HN = NewHNStructure();
  assert(HN);

  assert(CellMeshInterpretor::GetIntegrator()==NULL);
  BasicMeshInterpretor::GetIntegratorPointer() =  new GalerkinGlsIntegrator<3>;

  assert(CellMeshInterpretor::GetFem()==NULL);
  typedef Transformation3d<BaseQ13d>           TransQ1;
  typedef FiniteElement<3,2,TransQ1,BaseQ13d>  FiniteElement;
  BasicMeshInterpretor::GetFemPointer() =  new FiniteElement;

  CellMeshInterpretor::BasicInit(pf);
}

/* ----------------------------------------- */
