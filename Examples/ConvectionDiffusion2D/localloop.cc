#include  "localloop.h"
#include  "backup.h"
#include  "monitoring.h"
#include  "localmeshagent.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

void LocalLoop::BasicInit(const ParamFile* paramfile)
{
  LPD.BasicInit(paramfile);
  GetMeshAgentPointer() = new LocalMeshAgent;
  StdLoop::BasicInit(paramfile);
}

/*-------------------------------------------------*/

void LocalLoop::run()
{
  _iter=1;
  
  MultiLevelGhostVector u("u"), f("f");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  GlobalVector  ualt;

  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);

  Monitoring  Moning;

  vector<GlobalVector>  dat;

  MultiLevelSolverInterface* MP = GetMultiLevelSolver();

  _clock_newmesh.start();
  GetMultiLevelSolver()->ReInit(LPD);
  _clock_newmesh.stop();

  int nlevels = MP->nlevels();
  dat.resize(nlevels);
  for(int l=nlevels-1;l>=0;l--)
    {
      if(l==nlevels-1)
	{
	  string filename("../NavierStokes2D/Results/solve.00003.bup");
	  GlobalVector& d = dat[l];
	  ReadBackUpResize(d,filename);
	}
      else
	{
	  dat[l].ncomp() = dat[l+1].ncomp();
	  MP->GetSolver(l)->ResizeVector(&dat[l],"");
 	  MP->Transfer(l+1,dat[l],dat[l+1]);
	}
        
      MP->GetSolver(l)->AddNodeVector("beta",&dat[l]);
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
