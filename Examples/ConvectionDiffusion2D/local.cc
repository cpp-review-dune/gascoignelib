#include  "local.h"
#include  "backup.h"
#include  "monitoring.h"
#include  "compose_name.h"

using namespace std;
using namespace Gascoigne;

/*-------------------------------------------------*/

void LocalLoop::run(const ProblemDescriptorInterface* PD)
{
  _iter=1;
  
  VectorInterface u("u"), f("f"), dat("dat");
  GlobalVector  ualt;

  Monitoring  Moning;

  MultiLevelSolverInterface* MP = GetMultiLevelSolver();

  _clock_newmesh.start();
  GetMultiLevelSolver()->ReInit(*PD);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(f);
  GetMultiLevelSolver()->ReInitVector(dat);

  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  _clock_newmesh.stop();

  int nlevels = MP->nlevels();
  for(int l=nlevels-1;l>=0;l--)
    {
      GlobalVector& d = MP->GetSolver(l)->GetGV(dat);
      if(l==nlevels-1)
	{
	  string filename("solve.00003.bup");
	  ReadBackUpResize(d,filename);
	}
      else
	{
	  d.ncomp() = MP->GetSolver(l+1)->GetGV(dat).ncomp();
	  MP->GetSolver(l)->ReInitVector(dat);
 	  MP->Transfer(l+1,d,MP->GetSolver(l+1)->GetGV(dat));
	}
      MP->GetSolver(l)->AddNodeVector("beta",dat);
    }

  MP->GetSolver()->OutputSettings();
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
