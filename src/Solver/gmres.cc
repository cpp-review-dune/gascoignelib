#include  "../LinAlg/gmres.xx"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"

namespace Gascoigne
{
template GMRES<StdSolver,StdMultiLevelSolver,MultiLevelGhostVector>;
}
