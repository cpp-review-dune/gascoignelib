#include  "meshagent.h"
#include  "q1lps2d.h"
#include  "stdsolver.h"
#include  "stdmultilevelsolver.h"
#include  "onelevelalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "domainmeanfunctional.h"
#include  "problemdescriptorbase.h"
#include  "usefullfunctionsbd.h"
#include  "navierstokeslps2d.h"

using namespace Gascoigne;

/* ----------------------------------------- */

class BenchMarkDirichletData : public DirichletData
{
protected:
  double vmax;
public:
  BenchMarkDirichletData() {
    vmax = 0.3;
  }
  std::string GetName() const {return "Bench";}
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const {

    double y = v.y();

    b.zero();
    if (color!=80)
      {
        double high = 4.1;
        b[1] = vmax * ParabelFunction(y,0.,high);
      }
  }
};

/*---------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "NavierStokesBenchmark";}
  void BasicInit(const ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new NavierStokesLps2d(GetParamFile());
    GetDirichletDataPointer() = new BenchMarkDirichletData;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface*    NewDiscretization(int level) const { return new Q1Lps2d; }
  SolverInterface*            NewSolver(int level)         const { return new StdSolver;}
  MultiLevelSolverInterface*  NewMultiLevelSolver()        const { return new StdMultiLevelSolver;}
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
  PC.AddProblem("NavierStokesBenchmark", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  MultiLevelAlgorithm B;
  B.BasicInit(&paramfile,&N,&PC);
  B.GlobalRefineLoop("NavierStokesBenchmark");

  return 0;
}

