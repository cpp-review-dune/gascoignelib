#include  "stdloop.h"
#include  "starter.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"
#include  "laplace2d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "meshagent.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */
class LocalEquation : public Laplace2d
{
public:
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

  int dim=2;
  int prerefine=3;
  string inpname("square.inp");
  dynamic_cast<MeshAgent*>(GetMeshAgent())->BasicInit(dim, inpname, prerefine);

  StdLoop::BasicInit(paramfile);
  }
};

/*---------------------------------------------------*/
class LocalExactSolution : public ExactSolution
{
public:
  LocalExactSolution() : ExactSolution() {}

  std::string GetName() const {return "LocalExactSolution";}
  double operator()(int c, const Vertex2d& v)const{return v.x()*v.y()+11;}
  int GetNcomp() const { return 1; }
};

/*---------------------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorInterface
{
public:


private:


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
  
public:
  ProblemDescriptor() {}
  std::string GetName() const {return "Local";}
};

/*---------------------------------------------------*/

class LocalDomainFunctional : public virtual AllDomainFunctional
{
 public:

  LocalDomainFunctional() : AllDomainFunctional(1,0)
    {
      ExactValue() = 0.02776989201546093;
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
