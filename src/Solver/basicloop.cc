#include  "basicloop.h"
#include  "meshagent.h"
#include  "compose_name.h"
#include  "backup.h"
#include  "adaptordata.h"
#include  "filescanner.h"
#include  "monitoring.h"
#include  <iomanip>
#include  "gostream.h"
#include  "stdmultilevelsolver.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  BasicLoop::BasicLoop() : _paramfile(NULL), _MA(NULL), _ML(NULL), _iter(0), IOM("Results"), _SI(NULL)
{
  _reload  = "none";
}

/*-----------------------------------------*/

BasicLoop::~BasicLoop()
{
  if(_MA  != NULL) {delete _MA; _MA=NULL;}
  if(_ML  != NULL) {delete _ML; _ML=NULL;}
}

/*-----------------------------------------*/

void BasicLoop::ClockOutput() const
{
  cout << "************************************************************************\n\n";
  cout << "BasicLoop\t\tTIME\n";
  cout << "  NewMesh\t\t" << _clock_newmesh.read() << endl;
  cout << "  Solve\t\t\t" << _clock_solve.read() << endl;
  cout << "  Write\t\t\t" << _clock_write.read() << endl;
}

/*-----------------------------------------*/

void BasicLoop::BasicInit(const ParamFile* paramfile)
{
  _paramfile = paramfile;
  string s_copy_param_file;

  DataFormatHandler DFH;
  DFH.insert("niter",              &_niter,             1);
  DFH.insert("initial",            &_initial,           "boundary");
  DFH.insert("reload",             &_reload,            "none");
  DFH.insert("writebupgup",        &_writeBupGup,       true);
  DFH.insert("resultsdir",         &_s_resultsdir,      "Results");
  DFH.insert("copy_param_file",    &s_copy_param_file,  "no");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");

  Mon.init(_paramfile,1);
  Mon.set_directory(_s_resultsdir);

  // copy paramfile to Results directory?
  if(s_copy_param_file!="no"){ 
    string s_paramfile= paramfile->GetName();

    if(s_copy_param_file=="yes"){
      string cp_cmd = "cp "+s_paramfile+" "+_s_resultsdir; 
      //cout << "CP CMD:" << cp_cmd << "\n"; 
      system(cp_cmd.c_str()); 
    }else if(s_copy_param_file=="yes-with-date"){
      string cp_cmd = "cp "+s_paramfile+" "+_s_resultsdir + "/" + s_paramfile + "-`date +\"%Y.%m.%d-%H%M\"`" ; 
      //cout << "CP CMD:" << cp_cmd << "\n"; 
      system(cp_cmd.c_str()); 
    }else {
      cerr << __FILE__ << ":" << __LINE__ << ": bad copy_param_file parameter value="<< s_copy_param_file << endl;
      cerr << __FILE__ << ":" << __LINE__ << ": possible values: no, yes, yes-with-date"<< endl;
    }
  }


  assert((_reload=="none") || (_initial=="file"));

  if ((_reload!="none") && (_initial!="file"))
    {
      cerr << "Please, add 'initial file' to Block Loop" << endl;
      _initial = "file";
    }
  assert((_reload!="none") || (_initial!="file"));

  if(GetMeshAgentPointer()==NULL)
    {
      GetMeshAgentPointer() = new MeshAgent;
    }
  GetMeshAgent()->BasicInit(_paramfile);

  if(GetMultiLevelSolverPointer()==NULL)
    {
      GetMultiLevelSolverPointer() = new StdMultiLevelSolver;
    }
  assert(GetMultiLevelSolver());

  GetMultiLevelSolver()->BasicInit(GetMeshAgent(),_paramfile);
  GetMultiLevelSolver()->SetMonitorPtr(&Mon);

  if (GetSolverInfos()==NULL)
    {
      GetSolverInfosPointer() = new SolverInfos;
    }
  assert(GetSolverInfos());
  GetSolverInfos()->BasicInit(_paramfile);
}

/*-------------------------------------------------------*/

void BasicLoop::PrintMeshInformation(int outputlevel) const
{
  cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes();
  cout << " " << GetMeshAgent()->ncells() << endl;
  
  if(outputlevel)
    {
      for(int l=0;l<GetMeshAgent()->nlevels();l++)
        {
          const MeshInterface* M = GetMeshAgent()->GetMesh(l);
          cout << l << " [n,c] " << M->nnodes() << " " << M->ncells() << endl;
        }
    }
}

/*-------------------------------------------------------*/

void BasicLoop::Output(const MultiLevelGhostVector& u, string name) const
{
  GetMultiLevelSolver()->GetSolver()->Visu(name,u.finest(),_iter);
//   GetMultiLevelSolver()->GetSolver()->VisuGrid(name,_iter);
  if(_writeBupGup)
  {   
    WriteMeshAndSolution(name,u);
  }
}

/*-------------------------------------------------*/

void BasicLoop::WriteMeshAndSolution(const string& filename, const MultiLevelGhostVector& u) const
{
  string name;
  name = filename;
//   name = filename + "_value";
  compose_name(name,_iter);
  GetMultiLevelSolver()->GetSolver()->Write(u.finest(),name);
  cout << "[" << name << ".bup]";

//   name = filename + "_mesh";
//   compose_name(name,_iter);
  GetMeshAgent()->write_gup(name);
  cout << " [" << name << ".gup]" << endl;
}

/*-------------------------------------------------*/

void BasicLoop::WriteSolution(const MultiLevelGhostVector& u) const
{
  _clock_write.start();
  string filename = _s_resultsdir+"/solution";
  compose_name(filename,_iter);
  GetMultiLevelSolver()->GetSolver()->Write(u.finest(),filename);
  cout << "[" << filename << ".bup]" << endl;
  _clock_write.stop();
}

/*-------------------------------------------------*/

void BasicLoop::WriteMesh() const
{
  _clock_write.start();
  string filename = _s_resultsdir+"/mesh";
  compose_name(filename,_iter);
  GetMeshAgent()->write_gup(filename);
  cout << " [" << filename << ".gup]" << endl;
  _clock_write.stop();
}

/*-------------------------------------------------*/

void BasicLoop::InitSolution(MultiLevelGhostVector& u)
{
  u.zero();

  if      (_initial=="analytic") GetMultiLevelSolver()->GetSolver()->SolutionInit(u);
  else if (_initial=="file")     GetMultiLevelSolver()->GetSolver()->Read(u,_reload);
  else if (_initial=="boundary") GetMultiLevelSolver()->GetSolver()->BoundaryInit(u);
  else
    {
      cerr << "BasicLoop::InitSolution():\npossible values for init: \"analytic\", \"file\", \"boundary\"\n";
      assert(0);
    }
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  GetMultiLevelSolver()->GetSolver()->SubtractMean(u);
  GetMultiLevelSolver()->GetSolver()->Visu(_s_resultsdir+"/solve",u,0);
}

/*-------------------------------------------------*/

string BasicLoop::Solve(MultiLevelGhostVector& u, MultiLevelGhostVector& f, string name)
{
  _clock_solve.start();

  f.zero();
  GetMultiLevelSolver()->GetSolver()->Rhs(f.finest());

  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f.finest());
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u.finest());

  string status = GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
  _clock_solve.stop();

  _clock_write.start();
  Output(u,name);
  _clock_write.stop();

  return status;
}

/*-------------------------------------------------*/

void BasicLoop::ComputeGlobalErrors(const MultiLevelGhostVector& u)
{
  GetMultiLevelSolver()->GetSolver()->ComputeError(u,_GlobalErr);
  if (_GlobalErr.size()>0)
    {
      cout.precision(6);
      cout << "\nGlobalErrors l2,h1,l8 " << _GlobalErr << endl;
    }
}

/*-------------------------------------------------------*/

void BasicLoop::CopyVector(GlobalVector& dst, MultiLevelGhostVector& src)
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

void BasicLoop::CopyVector(MultiLevelGhostVector& dst, GlobalVector& src)
{
  int nn = src.n();
  int cc = src.ncomp();

  GetMultiLevelSolver()->GetSolver()->GetGV(dst).ncomp() = cc;
  GetMultiLevelSolver()->GetSolver()->GetGV(dst).resize(nn);
  GetMultiLevelSolver()->GetSolver()->GetGV(dst).equ(1.,src);
}

/*-------------------------------------------------*/

void BasicLoop::run(const ProblemDescriptorInterface* PD)
{
  MultiLevelGhostVector u("u"), f("f");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  GlobalVector  ualt;

  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  
  Monitoring Moning;
  
  for (_iter=1; _iter<=_niter; _iter++)
    {
      cout << "\n================== " << _iter << " ================";
      PrintMeshInformation();
      Moning.SetMeshInformation(_iter,GetMeshAgent()->nnodes(),GetMeshAgent()->ncells());
      
      _clock_newmesh.start();

      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      GetMultiLevelSolver()->ReInit(*PD);
      GetMultiLevelSolver()->InterpolateSolution(u,ualt);
      GetMultiLevelSolver()->GetSolver()->Visu(_s_resultsdir+"/interpolate",u,_iter);

      _clock_newmesh.stop();

      if (_iter==1) 
	{
	  GetMultiLevelSolver()->GetSolver()->OutputSettings();
	  InitSolution(u);
	}

      Solve(u,f);
      ComputeGlobalErrors(u);
      
      if (_iter<_niter) 
	{
	  CopyVector(ualt,u);

	  GetMeshAgent()->global_refine(1);
	}
     }
}
}

/*-------------------------------------------------*/
