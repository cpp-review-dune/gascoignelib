#include  "localmultilevelsolver.h"
#include  "localsolver.h"

using namespace std;

/* ----------------------------------------- */

LocalMultiLevelSolver::LocalMultiLevelSolver() : StdMultiLevelSolver()
{
}

/* ----------------------------------------- */

SolverInterface* LocalMultiLevelSolver::NewSolver(int solverlevel) 
{
  return new LocalSolver;
}
