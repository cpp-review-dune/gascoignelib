#include  "local.h"
#include  "backup.h"
#include  "monitoring.h"
#include  "compose_name.h"

using namespace std;
using namespace Gascoigne;

/*-------------------------------------------------*/

void LocalLoop::run(const std::string& problemlabel)
{
  _iter=1;

  VectorInterface u("u"), f("f"), dat("dat");
  GlobalVector  ualt;

  Monitoring  Moning;

  _clock_newmesh.start();
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(f);
  GetMultiLevelSolver()->ReInitVector(dat,3);

  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  _clock_newmesh.stop();

  GlobalVector& d = GetMultiLevelSolver()->GetSolver()->GetGV(dat);
  string filename("convection.bup");
  ReadBackUpResize(d,filename);

  GetMultiLevelSolver()->AddNodeVector("beta",dat);

  GetMultiLevelSolver()->GetSolver()->OutputSettings();
  InitSolution(u);
  Moning.BasicInit(GetExactValues());

  cout << "\n================== " << _iter << "================";
  cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes();
  cout << " " << GetMeshAgent()->ncells() << endl;
  Moning.SetMeshInformation(_iter,GetMeshAgent()->nnodes(),GetMeshAgent()->ncells());

  Solve(u,f);
  ComputeGlobalErrors(u);

  _clock_functionals.start();
  nvector<double> juh = Functionals(u,f);
  _clock_functionals.stop();

  nvector<double> eta;

  _clock_estimate.start();
  if (_estimator!="none")
  {
    double est = Estimator(eta,u,f);
    Moning.SetSolutionInformation(_iter,juh,est);
  }
  _clock_estimate.stop();
}
