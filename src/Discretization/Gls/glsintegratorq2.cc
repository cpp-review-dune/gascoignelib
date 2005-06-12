#include  "glsintegratorq2.h"


using namespace std;

namespace Gascoigne
{
/* ----------------------------------------- */

template<int DIM>
GlsIntegratorQ2<DIM>::GlsIntegratorQ2()
{
  if (DIM==2)
    GlsIntegrator<DIM>::IF = new QuadGauss9;
  else
    GlsIntegrator<DIM>::IF = new HexGauss27;
  assert(GlsIntegrator<DIM>::IF);
}

/*-----------------------------------------------------------*/

template class GlsIntegratorQ2<2>;
template class GlsIntegratorQ2<3>;
}
