#include  "localtimeloop.h"
#include  "stdtimesolver.h"
#include  "localmeshagent.h"
#include  "stdmultilevelsolver.h"
#include  "backup.h"
#include  "compose_name.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

void LocalTimeLoop::BasicInit(const ParamFile* paramfile)
{
  GetMeshAgentPointer() = new LocalMeshAgent;
  GetMeshAgent()->BasicInit(paramfile);
  StdTimeLoop::BasicInit(paramfile);
}

/* ----------------------------------------- */

void LocalTimeLoop::NewMesh(const ProblemDescriptorInterface* PD)
{
  GetMultiLevelSolver()->ReInit(*PD);
  
  cout << "\nLocalTimeLoop::NewMesh(): Mesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();
}

/* ----------------------------------------- */

void LocalTimeLoop::AddNodeVector(string filename)
{

  int nlevels = GetMultiLevelSolver()->nlevels();
  dat_node.resize(nlevels);
  cerr << "adding vector: " << filename << " " << nlevels << endl;
  for(int l=nlevels-1;l>=0;l--)
    {
      if(l==nlevels-1)
	{
	  GlobalVector& d = dat_node[l];
	  ReadBackUpResize(d,filename);
	}
      else
	{
	  dat_node[l].ncomp() = dat_node[l+1].ncomp();
	  GetMultiLevelSolver()->GetSolver(l)->ResizeVector(&dat_node[l],"");
	  GetMultiLevelSolver()->Transfer(l+1,dat_node[l],dat_node[l+1]);
	}
      GetMultiLevelSolver()->GetSolver(l)->AddNodeVector("U",&dat_node[l]);
    }
}

/* ----------------------------------------- */

void LocalTimeLoop::DeleteNodeVector()
{

  int nlevels = GetMultiLevelSolver()->nlevels();
  cerr << "deleting vector" << endl;
  for(int l=nlevels-1;l>=0;l--)
    {
      GetMultiLevelSolver()->GetSolver(l)->DeleteNodeVector("U");
    }
}

/* ----------------------------------------- */

void LocalTimeLoop::ReInit(const ProblemDescriptorInterface* PD)
{
  GetMultiLevelSolver()->ReInit(*PD);
}

/* ----------------------------------------- */

void LocalTimeLoop::init(string name, int iter, const ProblemDescriptorInterface* PD)
{
  MultiLevelGhostVector u("u"), f("f"), ualt("ualt");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  ualt.SetMultiLevelSolver(GetMultiLevelSolver());
  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  GetMultiLevelSolver()->RegisterVector(ualt);

  GetMultiLevelSolver()->ReInit(*PD);

  cerr << "§§§§§§§§§§§§§§§§§§§§§§\n";
  //NewMesh();

//   StdMultiLevelSolver* MS = dynamic_cast<StdMultiLevelSolver*>(GetMultiLevelSolver());
//   assert(MS);
//   MS->MemoryVector();
  cerr << "§§§§§§§§§§§§§§§§§§§§§§\n";
  
  nvector<double> eta;

  StdTimeSolver* TSS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TSS);
  
  TimeInfoBroadcast();

  // Anfangswerte
  InitSolution(u,f);
  
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  _iter=iter;
  Output(u,name);
}

/* ----------------------------------------- */

void LocalTimeLoop::backward(string iname, string name, int first, int last, const ProblemDescriptorInterface* PD)
{
  MultiLevelGhostVector u("u"), f("f"), ualt("ualt");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  ualt.SetMultiLevelSolver(GetMultiLevelSolver());
  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  GetMultiLevelSolver()->RegisterVector(ualt);
  
  nvector<double> eta;
  
  GetMultiLevelSolver()->ReInit(*PD);

  GetMultiLevelSolver()->GetSolver()->Read(u,iname);

  StdTimeSolver* TSS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TSS);

  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();
  //GetMultiLevelSolver()->GetSolver()->Visu("Results/solve",u,last);

  double T = 2000.;
  info.ReInitBackward(last,T);
  TimeInfoBroadcast();
  for (_iter=last; _iter>=first; _iter--)
    {
      string aname(name);
      compose_name(aname,_iter+1);
      aname += ".bup";
      cerr << aname << endl;
      AddNodeVector(aname);
			
      //
      // rhs fuer alten Zeitschritt
      //
      f.zero();
      GetMultiLevelSolver()->GetSolver()->TimeRhs(f,u);
      
      // neuer Zeitschritt
      //
      info.iteration_backward(_iter);
      TimeInfoBroadcast();

			ualt.equ(1.,u);

      SolveTimePrimal(u,f);
      Output(u,"Results/backward");
      Functionals(u,f);

      DeleteNodeVector();
    }
}

/* ----------------------------------------- */

void LocalTimeLoop::forward(string iname, int first, int last, const ProblemDescriptorInterface* PD)
{
  MultiLevelGhostVector u("u"), f("f"), ualt("ualt");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  ualt.SetMultiLevelSolver(GetMultiLevelSolver());
  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  GetMultiLevelSolver()->RegisterVector(ualt);
  
  nvector<double> eta;
  
  GetMultiLevelSolver()->ReInit(*PD);

  GetMultiLevelSolver()->GetSolver()->Read(u,iname);

  StdTimeSolver* TSS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TSS);

  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();
  //GetMultiLevelSolver()->GetSolver()->Visu("Results/solve",u,0);

  TimeInfoBroadcast();
  for (_iter=first; _iter<=last; _iter++)
    {
      // umschalten von Euler ?
      //
      info.SpecifyScheme(_iter);
			TimeInfoBroadcast();
      //
      // rhs fuer alten Zeitschritt
      //
      f.zero();
      GetMultiLevelSolver()->GetSolver()->TimeRhs(f,u);
      
      // neuer Zeitschritt
      //
      info.iteration(_iter);
      TimeInfoBroadcast();

      ualt.equ(1.,u);

      SolveTimePrimal(u,f
);
      Output(u,"Results/forward");
      Functionals(u,f);
    }
}
