#include "q2lps2dwithsecond.h"
#include "integratorlpswithsecond.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

void Q2Lps2dWithSecond::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    GetIntegratorPointer() =  new IntegratorLpsWithSecond<2>;
  assert(GetIntegrator());

  typedef FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond> FEWithSecond;

  if (GetFem()==NULL)
    GetFemPointer() =  new FEWithSecond;
  assert(GetFem());

  Q2Lps2d::BasicInit(paramfile);
}

/* ----------------------------------------- */
}
