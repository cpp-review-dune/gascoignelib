#include  "stdloop.h"
#include  "stdsolver.h"
#include  "stdmultilevelsolver.h"
#include  "q12d.h"
#include  "galerkinintegrator.h"
#include  "navierstokesgls2d.h"
#include  <string>
#include  "dirichletdata.h"
#include  "paramfile.h"
#include  "problemdescriptorinterface.h"
#include  "gascoignemeshconstructor.h"
#include  "levelmesh2d.h"
#include  "meshagent.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

class LocalDirichletData : public DirichletData
{
  protected:
    double vmax;
    
  public:
    LocalDirichletData(const ParamFile* paramfile)
    {
      vmax = 0.3;
    }
    ~LocalDirichletData() {}
    
    string GetName() const {return "Bench";}
    void operator()(Vector& b, const Vertex2d& v, int color) const 
    {
      b.zero();
      if (color!=4)
      {
        b[1] = vmax;
      }
    }
};

/*-----------------------------------------*/

class LocalEquation : public NavierStokesGls2d
{
  public:
    LocalEquation(const ParamFile* paramfile) : NavierStokesGls2d(paramfile) {}
    ~LocalEquation() {}
};

/*-----------------------------------------*/

class ProblemDescriptor : public ProblemDescriptorInterface
{
  protected:
  void ConstructEquation() 
  {
    GetEquationPointer() = new LocalEquation(GetParamFile());
  }
  void ConstructDirichletData() 
  {
    GetDirichletDataPointer() = new LocalDirichletData(GetParamFile());
  }
  
  public:
    ProblemDescriptor() : ProblemDescriptorInterface() {}
    ~ProblemDescriptor() {}

    string GetName() const {return "Local";}
};

/*-----------------------------------------*/

class LocalIntegrator : public GalerkinIntegrator<2>
{
  public:
    LocalIntegrator() : GalerkinIntegrator<2> () {}
    ~LocalIntegrator() {}
};

/*-----------------------------------------*/

class LocalMeshInterpretor : public Q12d
{
  public:
    LocalMeshInterpretor() : Q12d() {}
    ~LocalMeshInterpretor() {}
};

/*-----------------------------------------*/

class LocalSolver : public StdSolver
{
  public:
    LocalSolver() : StdSolver() {}
    ~LocalSolver() {}
    
//    MeshInterpretorInterface* NewMeshInterpretor(int dimension, const string& discname)
//    {
//     return new LocalMeshInterpretor;
//    }
};

/*-----------------------------------------*/

class LocalMultiLevelSolver : public StdMultiLevelSolver
{
  public:
    SolverInterface* NewSolver(int solverlevel)
    {
      return new LocalSolver;
    }
};

/*-----------------------------------------*/

/*-----------------------------------------*/
class LocalMeshAgent : public MeshAgent
{
  protected:
  public:
  LocalMeshAgent() : MeshAgent() {}
  ~LocalMeshAgent() {}
};

/*---------------------------------------------------*/

class LocalLoop : public StdLoop
{
  public:
    LocalLoop() : StdLoop() {}
    ~LocalLoop() {}
    
    void BasicInit(const ParamFile* paramfile) 
    {
      GetMeshAgentPointer() = new LocalMeshAgent;
      GetMeshAgent()->BasicInit(paramfile);
      GetMultiLevelSolverPointer() = new LocalMultiLevelSolver;
      StdLoop::BasicInit(paramfile);
    }
};
  
/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("bench.param");
  if(argc>=2) 
  {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  LocalLoop loop;
  loop.BasicInit(&paramfile);
  loop.run(&LPD);

  return 0;
}
