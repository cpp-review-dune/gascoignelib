#include  "localmultilevelsolver.h"
#include  "localsolver.h"


using namespace std;
using namespace Gascoigne;

/*-------------------------------------------------------------*/

SolverInterface* LocalMultiLevelSolver::NewSolver(int solverlevel) 
{
  return new LocalSolver;
}

/*-------------------------------------------------------------*/

void LocalMultiLevelSolver::SolutionTransfer(int l, GlobalVector& uc, GlobalVector& uf) const
{
  const LocalSolver* Sf = dynamic_cast<const LocalSolver*>(GetSolver(l+1));
  const LocalSolver* Sc = dynamic_cast<const LocalSolver*>(GetSolver(l));
  assert(Sf);

  int nc = Sc->n();
  uc.ReInit(uf.ncomp(),nc);

  Sf->HNAverage(uf);
  Interpolator[l]->SolutionTransfer(uc,uf);
  Sf->HNZero(uf);
//   Sc->SetBoundaryVector(uc,uc);
}
