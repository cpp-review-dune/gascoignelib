#include  "meshagent.h"
#include  "q12d.h"
#include  "stdsolver.h"
#include  "stdmultilevelsolver.h"
#include  "onelevelalgorithm.h"
#include  "multilevelalgorithm.h"
#include  "numericinterface.h"
#include  "domainmeanfunctional.h"
#include  "problemdescriptorbase.h"
#include  "laplace2d.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

class JumpDirichletData : public DirichletData
{
public:

  JumpDirichletData() {}
  std::string GetName() const {return "JumpDirichletData";}
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const 
  {
    double y = v.y();

    if (y<0.5) b[0] = 0.;
    else       b[0] = 1.;
  }
};

/*---------------------------------------------------*/

class NonlinearEllipticEquation : public Laplace2d
{
 private:

  double eps;

 public:

  NonlinearEllipticEquation() : Laplace2d(), eps(100.) {}

  std::string GetName()  const { return "NonlinearEllipticEquation";}
 
  void Form(VectorIterator b, const FemFunction& U, 
	    const TestFunction& N) const
  {
    Laplace2d::Form(b,U,N);
    b[0] += eps * (1.-U[0].m())*U[0].m()*N.m();
  }
  void Matrix(EntryMatrix& A, const FemFunction& U, 
	      const TestFunction& M, const TestFunction& N) const
  {
    Laplace2d::Matrix(A,U,M,N);
    A(0,0) -= eps * U[0].m()     *M.m()*N.m();
    A(0,0) += eps * (1.-U[0].m())*M.m()*N.m();
  }
};

/*---------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "egal";}
  void BasicInit(const ParamFile* pf) {
    GetEquationPointer()      = new NonlinearEllipticEquation;
    GetDirichletDataPointer() = new JumpDirichletData;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*----------------------------------------------------------------------------*/

class Numeric : public NumericInterface
{
public:

  DiscretizationInterface*    NewDiscretization(int level) const { return new Q12d; }
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
  PC.AddProblem("elliptic-nl", &LPD);

  //////////////
  // Discretization etc
  //////////////

  Numeric N;

  cout << "=================================" << endl;
  cout << "Algorithm with one level ilu solver:" << endl;
  cout << "=================================" << endl;

  OneLevelAlgorithm A;

  A.BasicInit(&paramfile,&N,&PC);
  A.RunNonLinear("elliptic-nl");

  cout << "=================================" << endl;
  cout << "Algorithm with Multilevel solver:" << endl;
  cout << "=================================" << endl;

  MultiLevelAlgorithm B;
  B.BasicInit(&paramfile,&N,&PC);
  B.RunNonLinear("elliptic-nl");

  return 0;
}

