#include  "newton.xx"
#include  "stdmultilevelsolver.h"
#include  "multilevelghostvector.h"

namespace Gascoigne
{
template void newton(StdMultiLevelSolver& S, MultiLevelGhostVector& u, const MultiLevelGhostVector& f, MultiLevelGhostVector& r, MultiLevelGhostVector& w, NLInfo& info);
}
