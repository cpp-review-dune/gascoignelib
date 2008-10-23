#include  "algorithm.h"
#include  "meshagent.h"
#include  "compose_name.h"
#include  "backup.h"
#include  "adaptordata.h"
#include  "filescanner.h"
#include  "monitoring.h"
#include  <iomanip>
#include  "gostream.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  Algorithm::Algorithm() : _MA(NULL), _SI(NULL), _paramfile(NULL)
{
}

/*-----------------------------------------*/

Algorithm::~Algorithm()
{
  if(_MA  != NULL) {delete _MA; _MA=NULL;}
  if(_SI  != NULL) {delete _SI; _SI=NULL;}
}

/*-----------------------------------------*/

void Algorithm::BasicInit(const ParamFile* paramfile, const NumericInterface* NI)
{
  _paramfile = paramfile;
  _NI = NI;

  GetMeshAgentPointer() = GetNumeric()->NewMeshAgent();
  GetMeshAgent()->BasicInit(_paramfile);

  GetSolverInfosPointer() = new SolverInfos;
  GetSolverInfos()->BasicInit(_paramfile);

  string command("mkdir -p Results");
  system(command.c_str());
}

/*-------------------------------------------------------*/

void Algorithm::PrintMeshInformation() const
{
  cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes();
  cout << " " << GetMeshAgent()->ncells() << endl;
}

/*-----------------------------------------*/

void Algorithm::Newton(VectorInterface& u, const VectorInterface& f, NLInfo& nlinfo)
{
  VectorInterface y("y"), du("du");

  ReInitVector(du); 
  ReInitVector(y); 

  cout << endl << "  ";
  nlinfo.reset();
  CGInfo& cginfo = nlinfo.GetLinearInfo();

  // Newton Residual

  GetSolver()->Equ(y,1.,f);
  GetSolver()->HNAverage(u);
  GetSolver()->Form(y,u,-1.);
  GetSolver()->HNZero(u);
  GetSolver()->SetBoundaryVectorZero(y);
  GetSolver()->SubtractMeanAlgebraic(y);

  double ny = GetSolver()->Norm(y);
  bool reached = nlinfo.check(0,ny,0.);

  for(int it=1; !reached; it++)
    {
      int nm1 = nlinfo.control().newmatrix();
      int nm2 = nlinfo.control().matrixmustbebuild();
      
      if (nm1+nm2!=0)
	{
	  if (nm1 && !nm2) cout << " N";
	  if (!nm1 && nm2) cout << "M ";
          if (nm1 && nm2)  cout << "MN";
	  AssembleMatrixAndIlu(u);
	  nlinfo.control().matrixmustbebuild() = 0;
	}
      else cout << "  ";
      GetSolver()->Zero(du);

      LinearSolve(du,y,cginfo);

      /////////////////////////////////
      //////// Newton update //////////

      double ndu = GetSolver()->Norm(du);

      if ( (cginfo.control().status()=="exploded") || 
	   (ndu>1.e10) || (!(ndu>=0.)) || (ny>1.e10) || (!(ny>=0.)) )
	{
	  nlinfo.control().status()="diverged";
	  cerr << "linear : " << cginfo.control().status() << endl;
	  cerr << "nonlinear : " << ndu << endl;
	}
      else
	{
	  double omega = 0.7;
	  double relax = 1.;
	  
	  for(int iter=0; iter<nlinfo.user().maxrelax(); iter++)
	    {
	      if(iter==0)
		{
		  GetSolver()->Add(u,relax,du);
		}
	      else
		{
		  GetSolver()->Add(u,relax*(omega-1.),du);
		  relax *= omega;
		}
	      GetSolver()->Equ(y,1.,f);
	      GetSolver()->HNAverage(u);
	      GetSolver()->Form(y,u,-1.);
	      GetSolver()->HNZero(u);
	      GetSolver()->SetBoundaryVectorZero(y);
	      GetSolver()->SubtractMeanAlgebraic(y);

	      ny = GetSolver()->Norm(y);

	      string message = nlinfo.check_damping(iter,ny);

	      if (message=="ok")       break;
	      if (message=="continue") continue;
	      if (message=="exploded") 
		{
		  GetSolver()->Add(u,-relax,du);
		  relax = 0.;
		  cout << "Damping exploded !!!!!" << endl;
		  nlinfo.control().status() = "diverged";
		  break;
		}
	    }
	}
      reached = nlinfo.check(it,ny,ndu);
    }
  DeleteVector(y);
  DeleteVector(du);
}

/*-------------------------------------------------------*/

void Algorithm::CopyVector(GlobalVector& dst, VectorInterface& src) const
{
  GetSolver()->HNAverage(src);
  
  int nn = GetSolver()->GetGV(src).n();
  int cc = GetSolver()->GetGV(src).ncomp();

  dst.ncomp() = cc;
  dst.resize(nn);
  
  dst.equ(1.,GetSolver()->GetGV(src));
  
  GetSolver()->HNZero(src);
}

/*-------------------------------------------------------*/
}
