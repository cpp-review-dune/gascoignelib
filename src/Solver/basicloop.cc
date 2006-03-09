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
  BasicLoop::BasicLoop() : _MA(NULL), _ML(NULL), _SI(NULL), _iter(0), _paramfile(NULL), IOM("Results") 
{
  _reload  = "none";
}

/*-----------------------------------------*/

BasicLoop::~BasicLoop()
{
  if(_MA  != NULL) {delete _MA; _MA=NULL;}
  if(_ML  != NULL) {delete _ML; _ML=NULL;}
  if(_SI  != NULL) {delete _SI; _SI=NULL;}
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

void BasicLoop::BasicInit(const ParamFile* paramfile,
			  const ProblemContainer* PC,
			  const FunctionalContainer* FC)
{
  _paramfile = paramfile;
  string s_copy_param_file;

  DataFormatHandler DFH;
  DFH.insert("niter",              &_niter,             1);
  DFH.insert("initial",            &_initial,           "boundary");
  DFH.insert("reload",             &_reload,            "none");
  DFH.insert("writevtk",           &_writeVtk,          true);
  DFH.insert("writebupgup",        &_writeBupGup,       true);
  DFH.insert("writeinp"   ,        &_writeInp   ,      false);  
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

  GetMultiLevelSolver()->BasicInit(GetMeshAgent(),_paramfile,PC,FC);
  GetMultiLevelSolver()->SetMonitorPtr(&Mon);

  if (GetSolverInfosPointer()==NULL)
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

void BasicLoop::Output(const VectorInterface& u, string name) const
{
  if(_writeVtk)
  {
    GetMultiLevelSolver()->GetSolver()->Visu(name,u,_iter);
  }
  if(_writeBupGup)
  {   
    WriteMeshAndSolution(name,u);
  }
  if(_writeInp)
  {   
    WriteMeshInp(name);
  }
}

/*-------------------------------------------------*/

void BasicLoop::WriteMeshAndSolution(const string& filename, const VectorInterface& u) const
{
  string name;
  name = filename;
//   name = filename + "_value";
  compose_name(name,_iter);
  GetMultiLevelSolver()->GetSolver()->Write(u,name);
  cout << "[" << name << ".bup]";

//   name = filename + "_mesh";
//   compose_name(name,_iter);
  GetMeshAgent()->write_gup(name);
  cout << " [" << name << ".gup]" << endl;
}

/*-------------------------------------------------*/

void BasicLoop::WriteSolution(const VectorInterface& u) const
{
  _clock_write.start();
  string filename = _s_resultsdir+"/solution";
  compose_name(filename,_iter);
  GetMultiLevelSolver()->GetSolver()->Write(u,filename);
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

void BasicLoop::WriteMeshInp(const string& name) const
{
  _clock_write.start();
  string filename = name;
  compose_name(filename,_iter);
  GetMeshAgent()->write_inp(filename);
  cout << " [" << filename << ".inp]" << endl;
  _clock_write.stop();
}

/*-------------------------------------------------*/

void BasicLoop::InitSolution(VectorInterface& u)
{
  GetMultiLevelSolver()->GetSolver()->Zero(u);

  if      (_initial=="analytic") GetMultiLevelSolver()->GetSolver()->SolutionInit(u);
  else if (_initial=="file")     GetMultiLevelSolver()->GetSolver()->Read(u,_reload);
  else if (_initial=="boundary") GetMultiLevelSolver()->GetSolver()->BoundaryInit(u);
  else if (_initial!="zero")
    {
      cerr << "BasicLoop::InitSolution():\npossible values for init: \"analytic\", \"file\", \"boundary\", \"zero\"\n";
      assert(0);
    }
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  GetMultiLevelSolver()->GetSolver()->SubtractMean(u);
  //  GetMultiLevelSolver()->GetSolver()->Visu(_s_resultsdir+"/solve",u,0);
}

/*-------------------------------------------------*/

string BasicLoop::Solve(VectorInterface& u, VectorInterface& f, string name)
{
  _clock_solve.start();

  GetMultiLevelSolver()->GetSolver()->Zero(f);
  GetMultiLevelSolver()->GetSolver()->Rhs(f);

  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);

  string status = GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
  _clock_solve.stop();

  _clock_write.start();
  Output(u,name);
  _clock_write.stop();

  return status;
}

/*-------------------------------------------------*/

void BasicLoop::ComputeGlobalErrors(const VectorInterface& u)
{
  GetMultiLevelSolver()->GetSolver()->ComputeError(u,_GlobalErr);
  if (_GlobalErr.size()>0)
    {
      cout.precision(6);
      for (int c=0; c<_GlobalErr.ncomp(); c++)
	{
	  cout << "\nGlobalErrors [" << c << "]  l2,h1,l8 ";
	  for (int i=0; i <_GlobalErr.n(); i++)
	    cout << " " << _GlobalErr(i,c);
	}
      cout << endl;
    }
}

/*-------------------------------------------------------*/

void BasicLoop::CopyVector(GlobalVector& dst, VectorInterface& src)
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

void BasicLoop::CopyVector(VectorInterface& dst, GlobalVector& src)
{
  int nn = src.n();
  int cc = src.ncomp();

  GetMultiLevelSolver()->GetSolver()->GetGV(dst).ncomp() = cc;
  GetMultiLevelSolver()->GetSolver()->GetGV(dst).resize(nn);
  GetMultiLevelSolver()->GetSolver()->GetGV(dst).equ(1.,src);
}

/*-------------------------------------------------*/

void BasicLoop::run(const std::string& problemlabel)
{
  VectorInterface u("u"), f("f");
  GlobalVector  ualt;

  Monitoring Moning;
  
  for (_iter=1; _iter<=_niter; _iter++)
    {
      cout << "\n================== " << _iter << " ================";
      PrintMeshInformation();
      Moning.SetMeshInformation(_iter,GetMeshAgent()->nnodes(),GetMeshAgent()->ncells());
      
      _clock_newmesh.start();

      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      GetMultiLevelSolver()->ReInit(problemlabel);
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(f);
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
