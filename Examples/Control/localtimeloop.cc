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
 	  GetMultiLevelSolver()->SolutionTransfer(l+1,dat_node[l],dat_node[l+1]);
	}
      GetMultiLevelSolver()->GetSolver(l)->AddNodeVector(&dat_node[l]);
    }
}

/* ----------------------------------------- */

void LocalTimeLoop::init(string name, int iter)
{
  NewMultiLevelGhostVector u("u"), f("f"), ualt("ualt");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  ualt.SetMultiLevelSolver(GetMultiLevelSolver());
  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  GetMultiLevelSolver()->RegisterVector(ualt);

  cerr << "§§§§§§§§§§§§§§§§§§§§§§\n";
  //NewMesh();

//   StdMultiLevelSolver* MS = dynamic_cast<StdMultiLevelSolver*>(GetMultiLevelSolver());
//   assert(MS);
//   MS->MemoryVector();
//   cerr << "§§§§§§§§§§§§§§§§§§§§§§\n";
  
  nvector<double> eta;

  StdTimeSolver* TSS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TSS);
  
  // Anfangswerte
  L2Projection(u,f);
  
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  _iter=iter;
  Output(u,name);
}

/* ----------------------------------------- */

void LocalTimeLoop::backward(string iname, string name, int first, int last, const ProblemDescriptorInterface* PD)
{
  NewMultiLevelGhostVector u("u"), f("f"), ualt("ualt");
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
  for (_iter=last; _iter>=first; _iter--)
    {
      info.iteration_backward(_iter);

      ualt.equ(1.,u);

      TimeInfoBroadcast();

      string aname(name);
      compose_name(aname,_iter);
      aname += ".bup";
      cerr << aname << endl;
      AddNodeVector(aname);

      SolveTimePrimal(u,f,"Results/backward");
      Functionals(u,f);
    }
}

/* ----------------------------------------- */

void LocalTimeLoop::forward(string iname, int first, int last, const ProblemDescriptorInterface* PD)
{
  NewMultiLevelGhostVector u("u"), f("f"), ualt("ualt");
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

  for (_iter=first; _iter<=last; _iter++)
    {
      info.iteration(_iter);

      ualt.equ(1.,u);

      TimeInfoBroadcast();

      SolveTimePrimal(u,f,"Results/forward");
      Functionals(u,f);
    }
}
