#include  "meshagent.h"
#include  "q1lps3d.h"
#include  "stdsolver.h"
#include  "stdmultilevelsolver.h"
#include  "onelevelalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "domainmeanfunctional.h"
#include  "problemdescriptorbase.h"
#include  "usefullfunctionsbd.h"
#include  "navierstokeslps3d.h"
#include  "boundaryfunction.h"

using namespace Gascoigne;

/* ----------------------------------------- */

class ThisDirichletData : public DirichletData
{
protected:

public:
  ThisDirichletData() 
  {
  }
  
  std::string GetName() const {return "Guten Tag";}
  void operator()(DoubleVector& b, const Vertex3d& v, int color) const 
  {
    double y = v.y();
    double z = v.z();

    b.zero();
    if (color==8)
      {
        b[1] = ParabelFunction(y,0.,1.0)*ParabelFunction(z,0.,1.0);
      }
  }
};

/*---------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Hallo";}
  void BasicInit(const ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new NavierStokesLps3d(GetParamFile());
    GetDirichletDataPointer() = new ThisDirichletData;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/* ----------------------------------------- */

class RundeKugel : public BoundaryFunction<3>
{
  double   _r;
  Vertex3d _c;
public :

  std::string GetName() const { return "Eine RundeKugel Eis";}
  void BasicInit(Vertex3d c, double r) {
      _c = c;
      _r = r;
    }
  double operator()(const Vertex3d& c) const {
      double r = - _r;
      for (int i=0; i<3; i++)
        {
          double dx = c[i]-_c[i];
          r += dx * dx;
        }
      return r;
    }
};

/*----------------------------------------------------------------------------*/

class CurvedMeshAgent : public MeshAgent
{
 protected:

  RundeKugel RK;

 public:

  CurvedMeshAgent() : MeshAgent()
    {
      double r = 0.05;
      Vertex3d v(0.5,0.5,0.5);
      RK.BasicInit(v,r);

      AddShape(80,&RK);
    }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface*    NewDiscretization(int level) const { return new Q1Lps3d; }
  SolverInterface*            NewSolver(int level)         const { return new StdSolver;}
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
  B.LocalRefineLoop ("NavierStokesBenchmark");

  return 0;
}

