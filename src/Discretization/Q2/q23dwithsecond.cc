#include "q23dwithsecond.h"
#include "integratorwithsecond.h"

using namespace std;
using namespace Gascoigne;

/**********************************************************/

void Q23dWithSecond::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    GetIntegratorPointer() =  new IntegratorWithSecond<3>;
  assert(GetIntegrator());

  typedef FiniteElementWithSecond<3, 2, Transformation3d<BaseQ23dWithSecond>, BaseQ23dWithSecond> FEWithSecond;

  if (GetFem()==NULL)
    GetFemPointer() =  new FEWithSecond;
  assert(GetFem());

  Q23d::BasicInit(paramfile);
}

/* ----------------------------------------- */

