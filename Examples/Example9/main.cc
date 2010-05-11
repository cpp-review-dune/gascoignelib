#include  "meshagent.h"
#include  "q22d.h"
#include  "timesolver.h"
#include  "multilevelsolver.h"
#include  "nonstationaryalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "domainmeanfunctional.h"
#include  "problemdescriptorbase.h"
#include  "boundaryfunction.h"
#include  "waveequation.h"
#include  "compose_name.h"
#include  "usefullfunctionsbd.h"
#include  "newmarkalgorithm.h"
#include  "newmarksolver.h"

using namespace Gascoigne;

/* ----------------------------------------- */

class WaveInitialCondition : public DomainInitialCondition
{
 public:
  
  WaveInitialCondition() : DomainInitialCondition() {}
  std::string GetName() const { return "WaveInitialCondition";}
  int        GetNcomp() const { return 1; }

  double operator()(int c, const Vertex2d& v)const 
  {
    //return StiffnessFunction(v.x(),0.5,20.);
    double r = (v.x()-0.5)*(v.x()-0.5)+(v.y()-0.5)*(v.y()-0.5);
    return exp(-50.*r);
  }
};

/*----------------------------------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "OceanCrustProblem";}
  void BasicInit(const ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new WaveEquation;
    GetInitialConditionPointer() = new WaveInitialCondition;
    GetBoundaryEquationPointer() = new WaveBoundaryEquation;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface* NewDiscretization(int level) const { return new Q22d; }
  SolverInterface*         NewSolver(int level)         const { return new NewmarkSolver;}
  MeshAgentInterface*      NewMeshAgent()               const { return new MeshAgent;}
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
  PC.AddProblem("Wave", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  MultiLevelSolver    S;
  NewmarkAlgorithm    B;
  B.BasicInit(&paramfile,&S,&N,&PC);
  B.Run("Wave");

  return 0;
}

