#include  "q2lps2d.h"
#include  "galerkinlpsintegratorq2.h"

namespace Gascoigne
{
/* ----------------------------------------- */

void Q2Lps2d::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    GetIntegratorPointer() =  new GalerkinLpsIntegratorQ2<2>;
  assert(GetIntegrator());

  Q22d::BasicInit(paramfile);
}

/* ----------------------------------------- */

}
