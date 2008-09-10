//#include  "heatproblem.h"
#include  "benchmarkproblem.h"
#include  "benchmarkmeshagent.h"
#include  "q1lps2d.h"
#include  "stdmultilevelsolver.h"
#include  "nonstationaryalgorithm.h"
#include  "numericinterface.h"
#include  "timesolver.h"

using namespace Gascoigne;

/*----------------------------------------------------------------------------*/

class TimeMultiLevelSolver : public StdMultiLevelSolver
{
public:
  
  SolverInterface* NewSolver(int solverlevel)   {  return new TimeSolver;}
  string GetName() const {  return "TimeMultiLevelSolver";}
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  //DiscretizationInterface*    NewDiscretization(int level) const { return new Q12d; }
  DiscretizationInterface*    NewDiscretization(int level) const { return new Q1Lps2d; }
  SolverInterface*            NewSolver(int level)         const { return new TimeSolver;}
  MultiLevelSolverInterface*  NewMultiLevelSolver()        const { return new TimeMultiLevelSolver;}
  MeshAgentInterface*         NewMeshAgent()               const { return new CurvedMeshAgent;}
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
  BenchMarkProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);
  ProblemContainer PC;
  PC.AddProblem("NavierStokesBenchmark", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  NonstationaryAlgorithm B;
  B.BasicInit(&paramfile,&N,&PC);
  //B.ImplicitEuler("NavierStokesBenchmark");
  //B.ThetaScheme("NavierStokesBenchmark");
  B.FractionalStepThetaScheme("NavierStokesBenchmark");

  return 0;
}

