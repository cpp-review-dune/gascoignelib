#include  "stdtimeloop.h"
#include  "filescanner.h"
#include  "stdtimesolver.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

void StdTimeLoop::BasicInit(const ParamFile* paramfile)
{
  StdLoop::BasicInit(paramfile);

  double tbegin, tend, deltat, theta;
  int    neuler;
  string scheme;

  DataFormatHandler DFH;
  DFH.insert("dt"    ,&deltat  ,1.);
  DFH.insert("tbegin",&tbegin  ,0.);
  DFH.insert("tend"  ,&tend    ,1.e4);
  DFH.insert("neuler",&neuler  ,10);
  DFH.insert("scheme",&scheme  ,"Euler");
  DFH.insert("theta" ,&theta   ,0.5);
  DFH.insert("reload",&_reload);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");

  _timeinfo.ReInit(tbegin,tend,deltat,scheme,neuler,theta);
}

/*-------------------------------------------------*/

string StdTimeLoop::SolveTimePrimal(MultiLevelGhostVector& u, MultiLevelGhostVector& f)
{
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);

  string status = GetMultiLevelSolver()->Solve(u,f);

  return status;
}

/*-------------------------------------------------*/

void StdTimeLoop::adaptive_run(const ProblemDescriptorInterface* PD)
{
  MultiLevelGhostVector u("u"), f("f");
  GlobalVector ualt;
  
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  
  nvector<double> eta;

  for (_iter=1; _iter<=_niter; _iter++)
    {
      GetMultiLevelSolver()->ReInit(*PD);
      TimeInfoBroadcast();

      GetMultiLevelSolver()->InterpolateSolution(u,ualt);

      if (_iter==1) 
        {
          GetMultiLevelSolver()->GetSolver()->OutputSettings();
          InitSolution(u,f);
        }

      cout << "\n================== " << _iter << "================";
      cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes();
      cout << " " << GetMeshAgent()->ncells() << endl;
      
      // umschalten von Euler ?
      //
      _timeinfo.SpecifyScheme(_iter);
      TimeInfoBroadcast();
      //
      // rhs fuer alten Zeitschritt
      //
      f.zero();
      GetMultiLevelSolver()->GetSolver()->TimeRhsOperator(f,u);
      GetMultiLevelSolver()->GetSolver()->TimeRhs(1,f);
      
      // neuer Zeitschritt
      //
      _timeinfo.iteration(_iter);
      TimeInfoBroadcast();

      GetMultiLevelSolver()->GetSolver()->TimeRhs(2,f);
      SolveTimePrimal(u,f);
      Output(u,"Results/solve");
      
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
      TS->SetTimeData(_timeinfo.dt(), _timeinfo.theta(), _timeinfo.time(), _timeinfo.oldrhs());
    }
}

/*-------------------------------------------------*/

void StdTimeLoop::InitSolution(MultiLevelGhostVector& u, MultiLevelGhostVector& f)
{
  if (_initial=="analytic") 
    {
      StdTimeSolver* TS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver());
      assert(TS);
      f.zero();
      TS->RhsL2Projection(f);
      TS->L2Projection(u,f);
      TS->Write(u.finest(),"Results/initialu");
    }
  else
    {
      BasicLoop::InitSolution(u);
    }
}

/*-------------------------------------------------*/

void StdTimeLoop::run(const ProblemDescriptorInterface* PD)
{
  MultiLevelGhostVector u("u"), f("f");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  
  nvector<double> eta;
  
  GetMultiLevelSolver()->ReInit(*PD);
  
  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();
  
  TimeInfoBroadcast();

  // Anfangswerte
  InitSolution(u,f);
  
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  //GetMultiLevelSolver()->GetSolver()->Visu("Results/solve",u.finest(),0);

  for (_iter=1; _iter<=_niter; _iter++)
    {
      // umschalten von Euler ?
      //
      _timeinfo.SpecifyScheme(_iter);
      TimeInfoBroadcast();
      //
      // rhs fuer alten Zeitschritt
      //
      f.zero();
      GetMultiLevelSolver()->GetSolver()->TimeRhsOperator(f,u);
      GetMultiLevelSolver()->GetSolver()->TimeRhs(1,f);
      
      // neuer Zeitschritt
      //
      _timeinfo.iteration(_iter);
      TimeInfoBroadcast();

      GetMultiLevelSolver()->GetSolver()->TimeRhs(2,f);

      SolveTimePrimal(u,f);
      Output(u,"Results/solve");

      Functionals(u,f);
    }
}
