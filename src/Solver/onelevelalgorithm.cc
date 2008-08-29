#include  "onelevelalgorithm.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

/*-----------------------------------------*/

void OneLevelAlgorithm::BasicInit(const ParamFile* paramfile, const NumericInterface* NI,
				  const ProblemContainer* PC)
{
  Algorithm::BasicInit(paramfile,NI);

  _S = GetNumeric()->NewSolver();

  GetSolver()->SetDiscretization(*GetNumeric()->NewDiscretization());
  GetSolver()->BasicInit(paramfile,GetMeshAgent()->GetDimension());

  _PC = PC;
}

/*-----------------------------------------*/

void OneLevelAlgorithm::IluSolver(VectorInterface& du, const VectorInterface& f, CGInfo& info)
{
  VectorInterface help("help");
  ReInitVector(help); 
  
  bool ok = info.check(GetSolver()->Norm(f),0.);
  for(int iter=0; !ok; iter++)
    {
      GetSolver()->MatrixResidual(help,du,f);
      double rnorm = GetSolver()->Norm(help);

      StdSolver* SS = dynamic_cast<StdSolver*>(GetSolver());
      SS->GetIlu()->solve(GetSolver()->GetGV(help));
      double cnorm = GetSolver()->Norm(help);

      GetSolver()->Add(du,1.,help);
      ok = info.check(rnorm,cnorm);
    }
  DeleteVector(help);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::RunLinear(const std::string& problemlabel)
{
  GetSolver()->NewMesh(GetMeshAgent()->GetMesh(0));
  
  GetSolver()->SetProblem(*_PC->GetProblem(problemlabel));

  GetSolver()->OutputSettings();
  PrintMeshInformation();

  VectorInterface u("u"), f("f"), du("du");

  ReInitVector(u); 
  ReInitVector(f); 
  ReInitVector(du); 

  GetSolver()->Zero(u); 
  GetSolver()->SetBoundaryVector(u);
  
  // Compute RHS
  
  GetSolver()->Zero(f); 
  GetSolver()->Rhs(f); 
  GetSolver()->SetBoundaryVector(f);
  
  // Compute Residual
  
  GetSolver()->HNAverage(u);
  GetSolver()->Form(f,u,-1.);
  GetSolver()->SetBoundaryVectorZero(f);
  
  // Assemble Matrix and ILU
  
  GetSolver()->RegisterMatrix();
  GetSolver()->ReInitMatrix();
  GetSolver()->MatrixZero();
  GetSolver()->AssembleMatrix(u,1.);
  GetSolver()->ComputeIlu();
  
  // Solve Linear System
  
  GetSolver()->Zero(du); 
  
  CGInfo& info = GetSolverInfos()->GetLInfo();
  
  IluSolver(du,f,info);

  cout  << endl << "Linear solver " << info.control().status() << endl << endl;
  GetSolver()->Add(u,1.,du);
  GetSolver()->HNZero(u);

  GetSolver()->Visu("Results/onelevel",u,0);

  DeleteVector(u);
  DeleteVector(f);
  DeleteVector(du);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::AssembleMatrixAndIlu(VectorInterface& u)
{
  GetSolver()->MatrixZero();
  GetSolver()->AssembleMatrix(u,1.);
  GetSolver()->ComputeIlu(u);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::LinearSolve(VectorInterface& du, const VectorInterface& y, CGInfo& cginfo)
{
  cginfo.reset();
  IluSolver(du,y,cginfo);
}

/*-----------------------------------------*/

void OneLevelAlgorithm::RunNonLinear(const std::string& problemlabel)
{
  GetSolver()->NewMesh(GetMeshAgent()->GetMesh(0));
  
  GetSolver()->SetProblem(*_PC->GetProblem(problemlabel));

  GetSolver()->OutputSettings();
  PrintMeshInformation();

  VectorInterface u("u"), f("f");

  ReInitVector(u); 
  ReInitVector(f); 

  GetSolver()->Zero(u); 
  GetSolver()->SolutionInit(u);;
  GetSolver()->SubtractMean(u);
  GetSolver()->SetBoundaryVector(u);
  
  // Compute RHS
  
  GetSolver()->Zero(f); 
  GetSolver()->Rhs(f); 
  GetSolver()->SetBoundaryVector(f);

  // Memory Matrix

  GetSolver()->RegisterMatrix();
  GetSolver()->ReInitMatrix();
  GetSolver()->MatrixZero();

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
  
  Newton(u,f,nlinfo);

  GetSolver()->Visu("Results/onelevel",u,0);

  DeleteVector(u);
  DeleteVector(f);
}
}
