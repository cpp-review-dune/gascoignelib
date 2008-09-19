#include  "benchsplitproblem.h"
#include  "../Example6/benchmarkmeshagent.h"
#include  "splittingsolver.h"
#include  "chorinalgorithm.h"
#include  "numericinterface.h"
#include  "timesolver.h"
#include  "q1lps2d.h"

using namespace Gascoigne;

/*----------------------------------------------------------------------------*/

class ChorinSolver : public SplittingSolver
{
protected:

  DiscretizationInterface* CreateDiscretization_1() const { return new Q1Lps2d;}
  DiscretizationInterface* CreateDiscretization_2() const { return new Q12d;}
};

/*----------------------------------------------------------------------------*/

class Numeric : public virtual NumericInterface
{
public:

  MeshAgentInterface*      NewMeshAgent()               const { return new CurvedMeshAgent;}
  SolverInterface*         NewSolver(int level)         const { return new ChorinSolver;}
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

  VelocityProblemDescriptor PV;
  PressureProblemDescriptor PP;
  PV.BasicInit(&paramfile);
  PP.BasicInit(&paramfile);
  ProblemContainer PC;
  PC.AddProblem("VelocityProblem", &PV);
  PC.AddProblem("PressureProblem", &PP);

  Numeric                N;
  MultiLevelSolver       S;
  ChorinAlgorithm        B;

  B.BasicInit(&paramfile,&S,&N,&PC);
  B.Run("VelocityProblem","PressureProblem");

  return 0;
}

