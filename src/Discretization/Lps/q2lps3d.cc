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
  if (GetIntegrator()==NULL)
    GetIntegratorPointer() =  new GalerkinLpsIntegratorQ2<3>;
  assert(GetIntegrator());

  Q23d::BasicInit(paramfile);
}
}
