#include  "glsintegratorq2.h"


using namespace std;

namespace Gascoigne
{
/* ----------------------------------------- */

template<int DIM>
GlsIntegratorQ2<DIM>::GlsIntegratorQ2()
{
  if (DIM==2)
    IF = new QuadGauss9;
  else
    IF = new HexGauss27;
  assert(IF);
}

/*-----------------------------------------------------------*/

template GlsIntegratorQ2<2>;
template GlsIntegratorQ2<3>;
}
