#include  "stdloop.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"
#include  "equation.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "meshagent.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */
class LocalEquation : public Equation
{
public:
  string GetName() const {return "Local";}
  int ncomp() const {return 1;}
  void OperatorStrong(Vector& b, const FemFunction& U) const {
    b[0] += U[0].m() - U[0].D();
  }
  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const {
    b[0] += U[0].m()*N.m() + U[0].x()*N.x() + U[0].y()*N.y();
  }
  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const {
    A(0,0) +=  M.m()*N.m() +  M.x()*N.x() + M.y()*N.y();
  }
};

/* ----------------------------------------- */
class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile) {
  GetMeshAgentPointer() = new MeshAgent;
  
  if(paramfile==NULL) 
    {
      int dim=2;
      int prerefine=3;
      string inpname("square.inp");
      dynamic_cast<MeshAgent*>(GetMeshAgent())->BasicInit(dim, inpname, prerefine);
    }
  else
    {
      dynamic_cast<MeshAgent*>(GetMeshAgent())->BasicInit(paramfile);
    }

  StdLoop::BasicInit(paramfile);
  }
};

/*---------------------------------------------------*/
class LocalExactSolution : public ExactSolution
{
public:
  LocalExactSolution() : ExactSolution() {}

  std::string GetName() const {return "LocalExactSolution";}
  double operator()(int c, const Vertex2d& v)const{return v.x()*v.y();}
//   double operator()(int c, const Vertex2d& v)const{return v.x()*v.y()+11.;}
  int GetNcomp() const { return 1; }
};

/*---------------------------------------------------*/
class ProblemDescriptor : public ProblemDescriptorInterface
{
protected:
  void ConstructEquation() {
    GetEquationPointer() = new LocalEquation;
  }
  void ConstructExactSolution() {
    GetExactSolutionPointer() = new LocalExactSolution();
  }
  void ConstructRightHandSideData() {
      GetRightHandSideDataPointer() = new RightHandSideDataByEquation(GetEquation(), GetExactSolution());
    }
  void ConstructDirichletData() {
    GetDirichletDataPointer() = new DirichletDataByExactSolution(GetExactSolution());
  }
  void ConstructBoundaryManager() {
    GetBoundaryManagerPointer() = new BoundaryManager(GetParamFile());
    GetBoundaryManager()->AddDirichlet(1,0);
    GetBoundaryManager()->AddDirichlet(2,0);
    GetBoundaryManager()->AddDirichlet(3,0);
    GetBoundaryManager()->AddDirichlet(4,0);
  }
 public:
  ProblemDescriptor() : ProblemDescriptorInterface() {}
  std::string GetName() const {return "Local";}
};

/*---------------------------------------------------*/
class LocalDomainFunctional : public virtual AllDomainFunctional
{
 public:

  LocalDomainFunctional() : AllDomainFunctional(1,0)
    {
      ExactValue() = 11.25;
      beautifulname = "LocalDomain";
    }
  ~LocalDomainFunctional() {}
};

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  const ParamFile* paramfile(NULL);

  if(argc==2)
    {
      paramfile = new ParamFile(argv[1]); 
    }
  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(paramfile);

  /////////////
  // Loop
  /////////////
  LocalLoop loop;
  loop.BasicInit(paramfile);

  /////////////
  // Functionals
  /////////////
  LocalDomainFunctional j1;
  std::vector<const Functional*> J(1);
  J[0] = &j1;
  loop.SetFunctionals(J);
  
  loop.run(&LPD);

  if(paramfile!=NULL)
    {
      delete paramfile;
    }
  return 0;
}
