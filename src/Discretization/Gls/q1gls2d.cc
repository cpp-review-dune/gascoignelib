#include  "q1gls2d.h"
#include  "galerkinglsintegrator.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq12d.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
void Q1Gls2d::BasicInit(const ParamFile* pf)
{
  assert(HN==NULL);
  HN = NewHNStructure();
  assert(HN);

  assert(CellDiscretization::GetIntegratorPointer()==NULL);
  CellDiscretization::GetIntegratorPointer() =  new GalerkinGlsIntegrator<2>;

  assert(CellDiscretization::GetFemPointer()==NULL);
  typedef Transformation2d<BaseQ12d>           TransQ1;
  typedef FiniteElement<2,1,TransQ1,BaseQ12d>  FiniteElement;
  CellDiscretization::GetFemPointer() =  new FiniteElement;

  CellDiscretization::BasicInit(pf);
}
}

/* ----------------------------------------- */
