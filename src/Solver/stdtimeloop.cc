#include  "stdtimeloop.h"
#include  "filescanner.h"
#include  "stdtimesolver.h"

using namespace std;

/*-----------------------------------------*/

void StdTimeLoop::BasicInit(const string& paramfile, const ProblemDescriptorInterface& PD)
{
  StdLoop::BasicInit(paramfile, PD);

  double tbegin, tend, deltat, theta;
  int    neuler;
  string scheme;

  DataFormatHandler DFH;
  DFH.insert("dt"    ,&deltat  ,1.);
  DFH.insert("tbegin",&tbegin  ,0.);
  DFH.insert("tend"  ,&tend    ,1.e4);
  DFH.insert("neuler"  ,&neuler    ,10);
  DFH.insert("scheme" ,&scheme   ,"Euler");
  DFH.insert("theta" ,&theta   ,0.5);
  DFH.insert("reload",&_reload);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");

  info.init(tbegin,tend,deltat,scheme,neuler,theta);
}

/*-------------------------------------------------*/

string StdTimeLoop::SolveTimePrimal(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f, string name)
{
  f.zero();
  GetMultiLevelSolver()->GetSolver()->TimeRhs(f,u);  // hier noch u=ualt !!

  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);

  string status = GetMultiLevelSolver()->Solve(u,f);

  GetMultiLevelSolver()->GetSolver()->Visu(name,u,_iter);
  Output(u,name);

  return status;
}

/*-------------------------------------------------*/

void StdTimeLoop::adaptive_run()
{
  NewMultiLevelGhostVector u("u"), f("f");
  CompVector<double> ualt;
  
  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  
  nvector<double> eta;
  
  for (_iter=1; _iter<=_niter; _iter++)
    {
      cout << "\n================== " << _iter << "================";
      cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes();
      cout << " " << GetMeshAgent()->ncells() << endl;
      
      info.iteration(_iter);
      GetMultiLevelSolver()->NewMesh(GetProblemDescriptor());
      
      GetMultiLevelSolver()->InterpolateSolution(u,ualt);

      if (_iter==1) 
	{
	  GetMultiLevelSolver()->GetSolver()->OutputSettings();
	  L2Projection(u,f);
	}
      TimeInfoBroadcast();

      SolveTimePrimal(u,f);
      
      GetMultiLevelSolver()->GetSolver()->EnergyEstimator(eta, u, f);

      cout << "eta " << eta.sum() << endl;

      if (_iter<_niter) 
	{
	  CopyVector(ualt,u);

	  AdaptMesh(eta);
	}
     }
}

/*-------------------------------------------------*/

void StdTimeLoop::TimeInfoBroadcast()
{
  for(int l=0;l<GetMultiLevelSolver()->nlevels();l++)
    {
      StdTimeSolver* TS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver(l));
      assert(TS);
      TS->SetTimeData(info.dt(), info.theta(), info.time());
    }
}

/*-------------------------------------------------*/

void StdTimeLoop::L2Projection(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f)
{
  StdTimeSolver* TS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TS);
  f.zero();
  TS->RhsL2Projection(f,u);
  TS->L2Projection(u,f);
}

/*-------------------------------------------------*/

void StdTimeLoop::run()
{
  NewMultiLevelGhostVector u("u"), f("f"), ualt("ualt");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  ualt.SetMultiLevelSolver(GetMultiLevelSolver());
  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  GetMultiLevelSolver()->RegisterVector(ualt);
  
  nvector<double> eta;
  
  GetMultiLevelSolver()->NewMesh(GetProblemDescriptor());

  StdTimeSolver* TSS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TSS);
  
  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();
  
  // Anfangswerte
  L2Projection(u,f);
  
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  GetMultiLevelSolver()->GetSolver()->Visu("Results/solve",u,0);

  for (_iter=1; _iter<_niter; _iter++)
    {
      info.iteration(_iter);

      ualt.equ(1.,u);

      TimeInfoBroadcast();

      SolveTimePrimal(u,f);
      Functionals(u,f);
    }
}
