#include  "q2lps2d.h"
#include  "galerkinlpsintegratorq2.h"

namespace Gascoigne
{
/* ----------------------------------------- */

void Q2Lps2d::BasicInit(const Gascoigne::ParamFile* paramfile)
{
  assert(GetIntegrator()==NULL);
  GetIntegratorPointer() =  new GalerkinLpsIntegratorQ2<2>;
  
  Q22d::BasicInit(paramfile);
}

/* ----------------------------------------- */

}
