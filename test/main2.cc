#include  "stdloop.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"
#include  "stokesgls2d.h"
#include  "navierstokesgls2d.h"
#include  "zerodirichletdata.h"
#include  "righthandsidedata.h"
#include  "meshagent.h"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"
#include  "q1gls2d.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */
class LocalEquation : public NavierStokesGls2d
{
public:
  string GetName() const {return "Local";}
  LocalEquation() : NavierStokesGls2d() {
  penalty = 0.; visc = 0.001;
  ST.delta0 = 1.;
  ST.alpha0 = 1.;
  ST.xeta0 = 6.; 
  cerr << "----------------------------------------------\n";
  cerr << visc << endl;
}

};

/* ----------------------------------------- */
class LocalNeumannData : public NeumannData
{
public:
  string GetName() const {return "Local";}
  int GetNcomp() const {return 3;}
  void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Vertex2d& v, const Vertex2d& n, int col) const {

    //b.zero();
  //assert(b.size()==3);

    if (col==1)
      {
	b[1] += -0.005 * n.x() * N.m();
	b[2] += -0.005 * n.y() * N.m();
	return;
      }
    if (col==2)
      {
	b[1] += -0.02 * n.x() * N.m();
	b[2] += -0.02 * n.y() * N.m();
	return;
      }
  }
};

/* ----------------------------------------- */
class LocalSolver : public StdSolver
{
  MeshInterpretorInterface* NewMeshInterpretor(int dimension, const string& discname) {
    return new Q1Gls2d;
  }
  
};
/* ----------------------------------------- */
class LocalMultiLevelSolver : public StdMultiLevelSolver
{
  SolverInterface* NewSolver(int solverlevel) {
    return new LocalSolver;
  }
};
/* ----------------------------------------- */
class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile) {
  GetMeshAgentPointer() = new MeshAgent;

  int dim=2;
  int prerefine=0;
  string inpname("mesh3.gup");
  dynamic_cast<MeshAgent*>(GetMeshAgent())->BasicInit(dim, inpname, prerefine);

  GetMultiLevelSolverPointer() = new LocalMultiLevelSolver;
  StdLoop::BasicInit(paramfile);
  }
};

/*---------------------------------------------------*/
class ProblemDescriptor : public ProblemDescriptorInterface
{
protected:
  void ConstructEquation() {
    GetEquationPointer() = new LocalEquation;
  }
  void ConstructDirichletData() {
    GetDirichletDataPointer() = new ZeroDirichletData();
  }
  void ConstructNeumannData() {
    GetNeumannDataPointer() = new LocalNeumannData();
  }
  void ConstructBoundaryManager() {
//     const ParamFile* paramfile = GetParamFile();
//     GetBoundaryManagerPointer() = new BoundaryManager(paramfile);

    const ParamFile* paramfile(NULL);
    GetBoundaryManagerPointer() = new BoundaryManager(paramfile);
    GetBoundaryManager()->AddDirichlet(3,1);
    GetBoundaryManager()->AddDirichlet(3,2);
    GetBoundaryManager()->AddDirichlet(8,1);
    GetBoundaryManager()->AddDirichlet(8,2);

    GetBoundaryManager()->AddNeumann(1);
    GetBoundaryManager()->AddNeumann(2);
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
