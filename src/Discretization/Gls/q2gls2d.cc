#include  "q2gls2d.h"
#include  "galerkinglsintegratorq2.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq12d.h"
#include  "baseq22d.h"
#include "integratorwithsecond.h"

using namespace std;
namespace Gascoigne
{
/* ----------------------------------------- */

void Q2Gls2d::BasicInit(const ParamFile* paramfile)
{
  assert(GetIntegrator()==NULL);
  GetIntegratorPointer() =  new GalerkinGlsIntegratorQ2<2>;
  GetIntegratorPointer()->BasicInit();

  assert(GetFem()==NULL);
  typedef Transformation2d<BaseQ12d>           TransQ1;
  typedef Transformation2d<BaseQ22d>           TransQ2;
  typedef FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond> FiniteElement;
  //typedef FiniteElement<2,1,TransQ2,BaseQ22d>  FiniteElement;
  PatchDiscretization::GetFemPointer() =  new FiniteElement;

  PatchDiscretization::BasicInit(paramfile);
}
}
