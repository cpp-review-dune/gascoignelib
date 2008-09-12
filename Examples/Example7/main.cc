#include  "../Example6/benchmarkproblem.h"
#include  "../Example6/benchmarkmeshagent.h"
#include  "q1lps2d.h"
#include  "stdmultilevelsolver.h"
#include  "dwralgorithm.h"
#include  "numericinterface.h"
#include  "timesolver.h"
#include  "benchmarkfunctionals.h"
#include  "functionalcontainer.h"
#include  "dirichletdatabycolor.h"

using namespace Gascoigne;

/*----------------------------------------------------------------------------*/

class DualDragProblemDescriptor : public ProblemDescriptorBase
{
 protected:

  const ResidualFunctional* FP;

 public:
  
  DualDragProblemDescriptor() : ProblemDescriptorBase(), FP(NULL) {}
  ~DualDragProblemDescriptor()
  {
    GetEquationPointer() = NULL;
    FP = NULL;
  }
  std::string GetName() const { return "DualDragProblemDescriptor";}

  void BasicInit(const ParamFile* filename, const Equation* equation, const ResidualFunctional* fp)
  {
    ProblemDescriptorBase::BasicInit(filename);
    
    GetEquationPointer() = const_cast<Equation*>(equation);

    FP = fp;

    int    comp  = FP->GetComps()[0];
    double scale = FP->GetScales()[0];
    std::set<int> col = FP->GetColors();
    
    GetDirichletDataPointer()  = new DirichletDataByColor(comp,col,scale);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface*    NewDiscretization(int level) const { return new Q1Lps2d; }
  SolverInterface*            NewSolver(int level)         const { return new TimeSolver;}
  MultiLevelSolverInterface*  NewMultiLevelSolver()        const { return new StdMultiLevelSolver;}
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
  BenchMarkProblemDescriptor P1;
  P1.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("NavierStokesBenchmark", &P1);

  /////////////
  // Functional
  /////////////
    
  DragFunctional j;   

  ////////////////////
  // Adjoint Problem
  ////////////////////
    
  DualDragProblemDescriptor P2;
  P2.BasicInit(&paramfile,P1.GetEquation(),&j);
  PC.AddProblem("DualDrag", &P2);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  DwrAlgorithm B;
  B.BasicInit(&paramfile,&N,&PC);
  B.AdaptiveLoop("NavierStokesBenchmark","DualDrag",j);

  return 0;
}

