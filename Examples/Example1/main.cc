#include  "meshagent.h"
#include  "q12d.h"
#include  "stdsolver.h"
#include  "stdmultilevelsolver.h"
#include  "onelevelalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "problemdescriptorbase.h"
#include  "laplace2d.h"
#include  "zerodirichletdata.h"
#include  "constantrighthandside.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const { return "egal";}

  void BasicInit(const ParamFile* pf) 
  {
    GetEquationPointer()      = new Laplace2d;
    GetRightHandSidePointer() = new ConstantRightHandSideData(1,0,123.);
    GetDirichletDataPointer() = new ZeroDirichletData;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface*    NewDiscretization(int level) const { return new Q12d; }
  SolverInterface*            NewSolver(int level)         const { return new StdSolver;}
  MeshAgentInterface*         NewMeshAgent()               const { return new MeshAgent;}
};

/*----------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("gascoigne.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);
  ProblemContainer PC;
  PC.AddProblem("laplace", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  cout << "=================================" << endl;
  cout << "Algorithm with one level ilu solver:" << endl;
  cout << "=================================" << endl;

  OneLevelAlgorithm A;

  A.BasicInit(&paramfile,&N,&PC);
  A.RunLinear("laplace");

  cout << "=================================" << endl;
  cout << "Algorithm with Multilevel solver:" << endl;
  cout << "=================================" << endl;

  MultiLevelSolver    S;
  MultiLevelAlgorithm B;
  B.BasicInit(&paramfile,&S,&N,&PC);
  B.RunLinear("laplace");

  return 0;
}

/*----------------------------------------------------------------------------*/
