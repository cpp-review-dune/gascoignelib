#include  "chorinalgorithm.h"
#include  "stdtimesolver.h"
#include  "compose_name.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

/*-------------------------------------------------*/

void ChorinAlgorithm::Run(const std::string& problem1, const std::string& problem2)
{
  DataFormatHandler DFH;
  DFH.insert("initial", &initial,"boundary");
  DFH.insert("dt"     , &dt    ,1.);
  DFH.insert("theta"  , &theta ,0.5);
  DFH.insert("niter"  , &niter ,1);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");

  GetMultiLevelSolver()->ReInit(problem1);
  GetMultiLevelSolver()->ReInit(problem2);

  GetSolver()->OutputSettings();

  TimeInfoBroadcast();

  if      (theta==1. ) ChorinUzawa(problem1,problem2);
  else if (theta==0.5) VanKan(problem1,problem2);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::ChorinUzawa(const std::string& problem1, const std::string& problem2)
{
  GetMultiLevelSolver()->SetProblem(problem1);
  VectorInterface v("v"), fv("fv"), u("u");
  ReInitVector(v);
  ReInitVector(u);
  ReInitVector(fv);
  InitSolution(initial,v);
  GetMultiLevelSolver()->AssembleMatrix(v);

  GetMultiLevelSolver()->SetProblem(problem2);
  VectorInterface p("p"), fp("fp"), pold("pold"), q("q");
  ReInitVector(p);
  ReInitVector(q);
  ReInitVector(pold);
  ReInitVector(fp);
  GetSolver()->AddNodeVector("velocity",v);
  GetSolver()->Zero(p);
  GetSolver()->Zero(pold);
  AssembleMatrixAndIlu(p);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
  CGInfo& cginfo = GetSolverInfos()->GetLInfo();

  double alpha = 0.5;
  double mu = 1.;

  for (int iter=1; iter<=niter; iter++)
    {
      GetMultiLevelSolver()->SetProblem(problem1);
      GetSolver()->Zero(fv);
      GetSolver()->AddNodeVector("pressure",p);
      GetSolver()->Rhs(fv,1.);
      GetSolver()->DeleteNodeVector("pressure");
      GetSolver()->AddNodeVector("pressure",pold);
      GetSolver()->Rhs(fv,-1.);
      GetSolver()->DeleteNodeVector("pressure");
      GetSolver()->MassMatrixVector(fv,v,1./dt);

      // neuer Zeitschritt
      //
      time += dt;
      TimeInfoBroadcast();

      cout << "\n============== " << iter << " ==== Chorin-Uzawa === ";
      cout << " [t,dt] "<< time << " " << dt << "\n";

      // velocity predictor
      //
      GetSolver()->SetBoundaryVector(v);
      GetSolver()->SetBoundaryVector(fv);
      nlinfo.reset();
      Newton(v,fv,nlinfo);
      //GetSolver()->Visu("Results/predictor",v,iter);

      // pressure Poisson problem
      //
      GetMultiLevelSolver()->SetProblem(problem2);
      GetSolver()->Zero(fp);
      GetSolver()->Zero(q);
      GetSolver()->Rhs(fp);
      LinearSolve(q,fp,cginfo);

      // velocity projection
      //
      GetMultiLevelSolver()->SetProblem(problem1);
      GetSolver()->Zero(fv);
      GetSolver()->MassMatrixVector(fv,v,1.);
      GetSolver()->AddNodeVector("pressure",q);
      GetSolver()->Rhs(fv,dt);
      GetSolver()->DeleteNodeVector("pressure");
      cginfo.reset();
      GetSplittingSolver()->InverseMassMatrix(v,fv,cginfo); 
      if (iter%20==0) GetSolver()->Visu("Results/velocity",v,iter);

      // pressure update
      //
      GetMultiLevelSolver()->SetProblem(problem2);
      GetSolver()->Equ(pold,1.,p);
      GetSolver()->Zero(fp);
      GetSolver()->Rhs(fp,alpha*mu*dt);
      GetSolver()->MassMatrixVector(fp,p,1.);
      cginfo.reset();
      GetSplittingSolver()->InverseMassMatrix(p,fp,cginfo); 
      if (iter%20==0) GetSolver()->Visu("Results/pressure",p,iter);
    }
  GetSolver()->DeleteNodeVector("velocity");
  DeleteVector(v);
  DeleteVector(u);
  DeleteVector(fv);
  DeleteVector(pold);
  DeleteVector(p);
  DeleteVector(q);
  DeleteVector(fp);
}

/*-------------------------------------------------*/

void ChorinAlgorithm::VanKan(const std::string& problem1, const std::string& problem2)
{
  GetMultiLevelSolver()->SetProblem(problem1);
  VectorInterface v("v"), fv("fv");
  ReInitVector(v);
  ReInitVector(fv);
  InitSolution(initial,v);
  GetMultiLevelSolver()->AssembleMatrix(v);

  GetMultiLevelSolver()->SetProblem(problem2);
  VectorInterface p("p"), fp("fp"), q("q");
  ReInitVector(p);
  ReInitVector(fp);
  ReInitVector(q);
  GetSolver()->AddNodeVector("velocity",v);
  GetSolver()->Zero(p);
  AssembleMatrixAndIlu(p);

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
  CGInfo& cginfo = GetSolverInfos()->GetLInfo();

  for (int iter=1; iter<=niter; iter++)
    {
      GetMultiLevelSolver()->SetProblem(problem1);
      GetSolver()->AddNodeVector("pressure",p);
      GetSolver()->Zero(fv);
      GetSolver()->Rhs(fv,1./theta);
      GetSolver()->DeleteNodeVector("pressure");

      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver());
      S->StdSolver::Form(fv,v,-(1.-theta)/theta);
      GetSolver()->MassMatrixVector(fv,v,1./(theta*dt));

      // neuer Zeitschritt
      //
      time += dt;
      TimeInfoBroadcast();

      cout << "\n============== " << iter << " ==== VanKan === ";
      cout << " [t,dt] "<< time << " " << dt << "\n";
      
      // velocity predictor
      //
      GetSolver()->SetBoundaryVector(v);
      GetSolver()->SetBoundaryVector(fv);
      nlinfo.reset();
      Newton(v,fv,nlinfo);
      GetSolver()->Visu("Results/predictor",v,iter);

      // pressure Poisson problem
      //
      GetMultiLevelSolver()->SetProblem(problem2);
      GetSolver()->Zero(fp);
      GetSolver()->Rhs(fp,1./theta);
      GetSolver()->Zero(q);
      LinearSolve(q,fp,cginfo);
      GetSolver()->Add(p,(1.-theta)/theta,q);
      GetSolver()->Visu("Results/pressure",p,iter);
//       string name = "Results/pressure";
//       compose_name(name,iter);
//       GetSolver()->Write(v,name);

      // velocity projection
      //
      GetMultiLevelSolver()->SetProblem(problem1);
      GetSolver()->Zero(fv);
      GetSolver()->AddNodeVector("pressure",q);
      GetSolver()->Rhs(fv,theta*dt);
      GetSolver()->DeleteNodeVector("pressure");
      GetSolver()->MassMatrixVector(fv,v,1.);
      GetSolver()->Zero(v);
      cginfo.reset();
      GetSplittingSolver()->InverseMassMatrix(v,fv,cginfo);      
      GetSolver()->Visu("Results/velocity",v,iter);
      
//       name = "Results/velocity";
//       compose_name(name,iter);
//       GetSolver()->Write(v,name);
    }
  GetSolver()->DeleteNodeVector("velocity");
  DeleteVector(v);
  DeleteVector(fv);
  DeleteVector(p);
  DeleteVector(fp);
}

}
