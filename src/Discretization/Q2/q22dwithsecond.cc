#include "q22dwithsecond.h"
#include "integratorwithsecond.h"

using namespace std;
using namespace Gascoigne;

/**********************************************************/

void Q22dWithSecond::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    GetIntegratorPointer() =  new IntegratorWithSecond<2>;
  assert(GetIntegrator());

  typedef FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond> FEWithSecond;

  if (GetFem()==NULL)
    GetFemPointer() =  new FEWithSecond;
  assert(GetFem());

  Q22d::BasicInit(paramfile);
}

/* ----------------------------------------- */

