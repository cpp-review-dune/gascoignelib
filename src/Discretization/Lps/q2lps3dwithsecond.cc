#include "q2lps3dwithsecond.h"
#include "integratorlpswithsecond.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

void Q2Lps3dWithSecond::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    GetIntegratorPointer() =  new IntegratorLpsWithSecond<3>;
  assert(GetIntegrator());

  typedef FiniteElementWithSecond<3, 2, Transformation3d<BaseQ23dWithSecond>, BaseQ23dWithSecond> FEWithSecond;

  if (GetFem()==NULL)
    GetFemPointer() =  new FEWithSecond;
  assert(GetFem());

  Q2Lps3d::BasicInit(paramfile);
}

/* ----------------------------------------- */
}
