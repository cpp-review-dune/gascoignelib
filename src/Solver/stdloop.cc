#include  "stdloop.h"
#include  "meshagent.h"
#include  "compose_name.h"
#include  "backup.h"
#include  "adaptordata.h"
#include  "diplomantenadaptor.h"
#include  "malteadaptor.h"
#include  "filescanner.h"
#include  "monitoring.h"
#include  <iomanip>
#include  "gostream.h"
#include  "stdmultilevelsolver.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

StdLoop::StdLoop() : _paramfile("none"), _MA(NULL), _ML(NULL), _FMP(NULL), _PD(NULL), _iter(0), IOM("Results")
{
  _estimator = _extrapolate = _reload  = "none";
}

/*-----------------------------------------*/

StdLoop::~StdLoop()
{
  cout << "************************************************************************\n\n";
  cout << "StdLoop\t\t\tTIME\n";
  cout << "  NewMesh\t\t" << _clock_newmesh.read() << endl;
  cout << "  Solve\t\t\t" << _clock_solve.read() << endl;
  cout << "  Functionals\t\t" << _clock_functionals.read() << endl;
  cout << "  Estimate\t\t" << _clock_estimate.read() << endl;
  cout << "  Write\t\t\t" << _clock_write.read() << endl;

  if(_MA  != NULL) {delete _MA; _MA=NULL;}
  if(_ML  != NULL) {delete _ML; _ML=NULL;}
  if(_FMP!= NULL) {delete _FMP; _FMP=NULL;}
}

/*-----------------------------------------*/

void StdLoop::BasicInit(const string& paramfile, const ProblemDescriptorInterface& PD)
{
  _PD = const_cast<ProblemDescriptorInterface*>(&PD);
  _paramfile = paramfile;

  DataFormatHandler DFH;
  DFH.insert("nmin",&_nmin,1000);
  DFH.insert("nmax",&_nmax,100000);
  DFH.insert("p",&_p,0.1);
  DFH.insert("coarse",&_coarse,0);
  DFH.insert("nmax",&_nmax,100000);
  DFH.insert("niter",&_niter,6);
  DFH.insert("refiner",&_refiner,"eta");
  DFH.insert("estimator",&_estimator,"energy");
  DFH.insert("extrapolate",&_extrapolate,"no");
  DFH.insert("initial",&_initial,"boundary");
  DFH.insert("reload",&_reload,"none");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");

  Mon.init(_paramfile,1);
  Mon.set_directory("Results");

  assert((_reload=="none") || (_initial=="file"));

  if ((_reload!="none") && (_initial!="file"))
    {
      cerr << "Please, add 'initial file' to Block StdLoop" << endl;
      _initial = "file";
    }
  assert((_reload!="none") || (_initial!="file"));

  if(GetMeshAgentPointer()==NULL)
    {
      GetMeshAgentPointer() = new MeshAgent;
    }
  assert(GetMeshAgent());
  GetMeshAgent()->BasicInit(_paramfile);

  if(GetMultiLevelSolverPointer()==NULL)
    {
      GetMultiLevelSolverPointer() = new StdMultiLevelSolver;
    }
  assert(GetMultiLevelSolver());

  //
  StdLoop::NewFunctionalManager();
  GetFunctionalManager()->ConstructSet(_paramfile,*GetProblemDescriptor()->GetEquation());
  //

  GetMultiLevelSolver()->BasicInit(GetMeshAgent(),_paramfile);
  GetMultiLevelSolver()->SetMonitorPtr(&Mon);
}

/*-------------------------------------------------------*/

void StdLoop::NewFunctionalManager()
{
  GetFunctionalManagerPointer() = new FunctionalManager;
  assert(GetFunctionalManager());
}

/*-------------------------------------------------------*/

nvector<double> StdLoop::ComputeAllFunctionals(NewMultiLevelGhostVector& f, NewMultiLevelGhostVector& u) const
{
  vector<string> fnames = GetFunctionalManager()->GetFunctionalNames();
  int n = fnames.size(); 
  nvector<double> j(n,0.);
  for(int i=0;i<n;i++)
    {
      const Functional* FP = GetFunctionalManager()->GetFunctional(fnames[i]);
      assert(FP);
      j[i] = GetMultiLevelSolver()->ComputeFunctional(f,u,FP);
    }
  return j;
} 

/*-----------------------------------------*/

void StdLoop::EtaVisu(string name, int i, const nvector<double>& eta) 
{
  Visualization Visu;
  Visu.format("vtk");
  Visu.set_name(name);
  Visu.step(i);
  VisuDataInfo        VDI(1);
  VisuDataNVector  VD(eta);
  
  Visu.SetPointData(&VD);
  Visu.SetMesh(*GetMeshAgent()->GetMesh(0));
  Visu.SetPointDataInfo(&VDI);
  Visu.write();
}

/*-------------------------------------------------*/

void StdLoop::Output(const NewMultiLevelGhostVector& u, string name) const
{
  GetMultiLevelSolver()->GetSolver()->Visu(name,u.finest(),_iter);
//   GetMultiLevelSolver()->GetSolver()->VisuGrid(name,_iter);
  WriteMeshAndSolution(name,u);
}

/*-------------------------------------------------*/

void StdLoop::WriteMeshAndSolution(const string& filename, const NewMultiLevelGhostVector& u) const
{
  string name;
  name = filename;
//   name = filename + "_value";
  compose_name(name,_iter);
  GetMultiLevelSolver()->GetSolver()->Write(u.finest(),name);
  cout << "[" << name << "]";

//   name = filename + "_mesh";
//   compose_name(name,_iter);
  GetMeshAgent()->write_gup(name);
  cout << " [" << name << "]" << endl;
}

/*-------------------------------------------------*/

void StdLoop::WriteSolution(const NewMultiLevelGhostVector& u) const
{
  _clock_write.start();
  string filename = "Results/solution";
  compose_name(filename,_iter);
  GetMultiLevelSolver()->GetSolver()->Write(u.finest(),filename);
  cout << "[" << filename << "]";
  _clock_write.stop();
}

/*-------------------------------------------------*/

void StdLoop::WriteMesh() const
{
  _clock_write.start();
  string filename = "Results/mesh";
  compose_name(filename,_iter);
  GetMeshAgent()->write_gup(filename);
  cout << " [" << filename << "]" << endl;
  _clock_write.stop();
}

/*-------------------------------------------------*/

void StdLoop::InitSolution(NewMultiLevelGhostVector& u)
{
  u.zero();

  if      (_initial=="analytic") GetMultiLevelSolver()->GetSolver()->SolutionInit(u);
  else if (_initial=="file")     GetMultiLevelSolver()->GetSolver()->Read(u,_reload);
  else if (_initial=="boundary") GetMultiLevelSolver()->GetSolver()->BoundaryInit(u);
  else
    {
      assert(0);
    }
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  GetMultiLevelSolver()->GetSolver()->Visu("Results/solve",u,0);
}

/*-------------------------------------------------*/

string StdLoop::Solve(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f, string name)
{
  _clock_solve.start();

  f.zero();
  GetMultiLevelSolver()->GetSolver()->Rhs(f.finest());

  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f.finest());
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u.finest());

  string status = GetMultiLevelSolver()->Solve(u,f);
  _clock_solve.stop();

  _clock_write.start();
  Output(u,name);
  _clock_write.stop();

  return status;
}

/*-------------------------------------------------------*/

nvector<double> StdLoop::GetExactValues() const
{
  vector<string> fnames = GetFunctionalManager()->GetFunctionalNames();
  int n = fnames.size(); 
  nvector<double> j(n,0.);
  for(int i=0;i<n;i++)
    {
      const Functional* FP = GetFunctionalManager()->GetFunctional(fnames[i]);
      assert(FP);
      j[i] = FP->ExactValue();
    }
  return j; 
}

/*-------------------------------------------------*/

nvector<double> StdLoop::Functionals(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f)
{
  nvector<double> J  = ComputeAllFunctionals(f,u);
  nvector<double> JE = GetExactValues();
  _JErr.resize(J.size());
  if (J.size())
    {
      cout << "\nFunctionals  ";
      cout.precision(16);
      for(int i=0; i<J.size(); i++) 
	{
	  cout << "\t" << J[i];
	}
      cout << "\nErrors  ";
      for(int i=0; i<J.size(); i++) 
	{
	  _JErr[i] = JE[i] - J[i];
	  cout << "\t" << _JErr[i] ;
	}
      cout << endl;
      if(_extrapolate=="yes")
	{
	  Extra.NewValues(J);
	  Extra.Print();
	}
      cout << endl;
    }
  return J;
}

/*-------------------------------------------------*/

double StdLoop::Estimator(nvector<double>& eta, NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f)
{
  vector<string> name = GetGridFunctionalNames();
  cout << "Estimator\t"; 

  cout << "Baustelle!!!!!\n";
  return 0.;

  double est = 0.;
  if (_estimator=="energy")
    {
      est = GetMultiLevelSolver()->GetSolver()->EnergyEstimator(eta, u, f);
    }
  else assert(0);

  double eff = 0.;
  if      ((_estimator=="weighted") && (_JErr.size()>0))      eff = est/_JErr[0];
  else if ((_estimator=="energy")   && (_GlobalErr.size()>0)) eff = est/_GlobalErr(1,0);
  if  (eff!=0.) cout << " @ " << eff << endl; 
  
  cout << endl; 
  EtaVisu("Results/eta",_iter,eta);
  return est;
}

/*-------------------------------------------------*/

void StdLoop::ComputeGlobalErrors(const NewMultiLevelGhostVector& u)
{
  GetMultiLevelSolver()->GetSolver()->ComputeError(u,_GlobalErr);
  cout << "\nGlobalErrors l2,h1,l8 " << _GlobalErr << endl;
}

/*-------------------------------------------------*/

void StdLoop::AdaptMesh(const nvector<double>& eta)
{
  if     (_refiner=="global") GetMeshAgent()->global_refine(1);
  else if(_refiner=="random") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
      GetMeshAgent()->random_patch_refine(_p,0);
    }
  else if(_refiner=="eta") 
    {
      nvector<int> refnodes, coarsenodes;

      MalteAdaptor A(_paramfile,eta);
      A.refine(refnodes,coarsenodes);

      GetMeshAgent()->refine_nodes(refnodes,coarsenodes);
    }
  else if(_refiner=="dip") 
    {
      nvector<int> refnodes, coarsenodes;

      AdaptorData info;
      info.rfactor() = 1.; 
      DiplomantenAdaptor A(info,eta);
      A.refine(refnodes);
      GetMeshAgent()->refine_nodes(refnodes,coarsenodes);
    }
  else assert(0);
}

/*-------------------------------------------------*/

void StdLoop::CopyVector(GlobalVector& dst, NewMultiLevelGhostVector& src)
{
  GetMultiLevelSolver()->GetSolver()->HNAverage(src);
  
  int nn = GetMultiLevelSolver()->GetSolver()->GetGV(src).n();
  int cc = GetMultiLevelSolver()->GetSolver()->GetGV(src).ncomp();

  dst.ncomp() = cc;
  dst.resize(nn);
  
  dst.equ(1.,GetMultiLevelSolver()->GetSolver()->GetGV(src));
  
  GetMultiLevelSolver()->GetSolver()->HNZero(src);
}

/*-------------------------------------------------*/

void StdLoop::CopyVector(NewMultiLevelGhostVector& dst, GlobalVector& src)
{
  int nn = src.n();
  int cc = src.ncomp();

  GetMultiLevelSolver()->GetSolver()->GetGV(dst).ncomp() = cc;
  GetMultiLevelSolver()->GetSolver()->GetGV(dst).resize(nn);
  GetMultiLevelSolver()->GetSolver()->GetGV(dst).equ(1.,src);
}

/*-------------------------------------------------*/

void StdLoop::run()
{
  NewMultiLevelGhostVector u("u"), f("f");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  GlobalVector  ualt;

  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  
  Monitoring Moning;
  
  for (_iter=1; _iter<=_niter; _iter++)
    {
      cout << "\n================== " << _iter << " ================";
      cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes();
      cout << " " << GetMeshAgent()->ncells() << endl;
      Moning.SetMeshInformation(_iter,GetMeshAgent()->nnodes(),GetMeshAgent()->ncells());
      
      _clock_newmesh.start();
      GetMultiLevelSolver()->NewMesh(GetProblemDescriptor());
      GetMultiLevelSolver()->InterpolateSolution(u,ualt);
      GetMultiLevelSolver()->GetSolver()->Visu("Results/interpolate",u,_iter);

      _clock_newmesh.stop();

      if (_iter==1) 
	{
	  GetMultiLevelSolver()->GetSolver()->OutputSettings();
	  InitSolution(u);
	  Moning.BasicInit(GetExactValues());
	}

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
      if (_iter<_niter) 
	{
	  CopyVector(ualt,u);

	  AdaptMesh(eta);
	}
      _clock_estimate.stop();
     }
}

/*-------------------------------------------------*/
