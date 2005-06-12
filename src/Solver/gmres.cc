#include  "../LinAlg/gmres.xx"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"

namespace Gascoigne
{
  template class GMRES<StdSolver,StdMultiLevelSolver,VectorInterface>;
}
