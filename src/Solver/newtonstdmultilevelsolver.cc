#include  "newton.xx"
#include  "stdmultilevelsolver.h"
#include  "newmultilevelghostvector.h"

template void newnewton(StdMultiLevelSolver& S, NewMultiLevelGhostVector& u, const NewMultiLevelGhostVector& f, NewMultiLevelGhostVector& r, NewMultiLevelGhostVector& w, NLInfo& info);
