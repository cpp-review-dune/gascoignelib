#include  "nonstationaryalgorithm.h"
#include  "stdtimesolver.h"
#include  "compose_name.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
/*-----------------------------------------*/

void NonstationaryAlgorithm::BasicInit(const ParamFile* paramfile, 
				       const NumericInterface* NI,
				       const ProblemContainer* PC)
{
  MultiLevelAlgorithm::BasicInit(paramfile,NI,PC);

  DataFormatHandler DFH;
  DFH.insert("dt"    ,&dt    ,1.);
  DFH.insert("theta" ,&theta ,1.);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"Loop");
  
  time = 0.;
}

/*-------------------------------------------------*/

void NonstationaryAlgorithm::TimeInfoBroadcast()
{
  for (int l=0; l<GetMultiLevelSolver()->nlevels(); l++)
    {
      GetSolver(l)->SetTimeData(dt,theta,time);
    }
}

/*-------------------------------------------------*/

void NonstationaryAlgorithm::InitSolution(const string& initial, VectorInterface& u) const
{
  GetMultiLevelSolver()->GetSolver()->Zero(u);

  std::string reload = "stationary.00003.bup";

  if      (initial=="analytic") GetSolver()->SolutionInit(u);
  else if (initial=="file")     GetSolver()->Read(u,reload);
  else if (initial=="boundary") GetSolver()->BoundaryInit(u);
  else if (initial!="zero")
    {
      cerr << "ERROR: NonstationaryAlgorithm::InitSolution" << endl;
      abort();
    }
  GetSolver()->SetBoundaryVector(u);
  GetSolver()->SubtractMean(u);
  GetSolver()->Visu("Results/solve",u,0);
}

/*-----------------------------------------*/

void NonstationaryAlgorithm::ImplicitEuler(const std::string& problemlabel)
{
  theta = 1.;
  ThetaScheme(problemlabel);
}

/*-----------------------------------------*/

void NonstationaryAlgorithm::CrankNicholson(const std::string& problemlabel)
{
  theta = 0.5;
  ThetaScheme(problemlabel);
}

/*-----------------------------------------*/

void NonstationaryAlgorithm::ThetaScheme(const std::string& problemlabel)
{
  int    niter;
  string initial;

  DataFormatHandler DFH;
  DFH.insert("niter", &niter, 1);
  DFH.insert("initial", &initial,"boundary");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");
  
  assert(theta>0.);

  GetMultiLevelSolver()->ReInit(problemlabel);
  GetSolver()->OutputSettings();

  TimeInfoBroadcast();

  VectorInterface u("u"), f("f");
  ReInitVector(u);
  ReInitVector(f);
  InitSolution(initial,u);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();

  time += dt;

  for (int iter=1; iter<=niter; iter++)
    {
      cout << "\n============== " << iter << " ==== theta-scheme === ";
      cout << " [t,dt] "<< time << " " << dt << "\n";

      //
      // rhs fuer alten Zeitschritt
      //
      GetSolver()->Zero(f);
      GetSolver()->Rhs(f);
      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver());
      if (theta!=1.) S->StdSolver::Form(f,u,1.-1./theta);
      GetSolver()->MassMatrixVector(f,u,1./(theta*dt));
      
      // neuer Zeitschritt
      //
      time += dt;
      TimeInfoBroadcast();

      if (theta!=1.) GetSolver()->Rhs(f,1./theta-1.);
      GetSolver()->SetBoundaryVector(f);

      nlinfo.reset();
  
      Newton(u,f,nlinfo);

      GetSolver()->Visu("Results/solve",u,iter);
    }
  DeleteVector(u);
  DeleteVector(f);
}

/*-----------------------------------------*/


}
