#include  "newton.xx"
#include  "stdmultilevelsolver.h"
#include  "multilevelghostvector.h"

template void newton(StdMultiLevelSolver& S, MultiLevelGhostVector& u, const MultiLevelGhostVector& f, MultiLevelGhostVector& r, MultiLevelGhostVector& w, NLInfo& info);
