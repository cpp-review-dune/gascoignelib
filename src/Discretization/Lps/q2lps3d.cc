#include  "q2lps3d.h"
#include  "galerkinlpsintegratorq2.h"
#include  "transformation3d.h"
#include  "finiteelement.h"
#include  "baseq23d.h"

namespace Gascoigne
{
/* ----------------------------------------- */

void Q2Lps3d::BasicInit(const ParamFile* paramfile)
{
  assert(GetIntegrator()==NULL);
  GetIntegratorPointer() =  new GalerkinLpsIntegratorQ2<3>;

  assert(GetFem()==NULL);
  typedef Transformation3d<BaseQ23d>           TransQ2;
  typedef FiniteElement<3,2,TransQ2,BaseQ23d>  FiniteElement;
  PatchDiscretization::GetFemPointer() =  new FiniteElement;

  PatchDiscretization::BasicInit(paramfile);
}
}
