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

class LocalMesh2d : public GascoigneMesh2d
{
  protected:
    int _count;
    vector<nvector<double> > _q;
    map<int,BoundaryFunction<2>*>  MyShapes;
    void boundary_newton2d(){}

  public:
    LocalMesh2d() 
    {
      _count=0;
      _q.resize(10);
      for(int k=0;k<_q.size();k++)
      {
        _q[k].reservesize(10);
        _q[k].zero();
        _q[k][k]=0.25;
//        for(int n=0;n<_q[k].size();n++)
//        {
//          _q[k][n] = 0.25*pow(1./(n+1.),k);
//        }
      }
    }
    ~LocalMesh2d() {}

    void SetCoordinates(const LevelMesh2d* LM) 
    {
      assert(_count<_q.size());
      
      double rect_min = 0.75;
      double rect_max = 3.;
      for(int i=0;i<nnodes();i++)
      {
        const Vertex2d& v = LM->vertex2d(i);
        double x = v.x();
        double y = v.y();
        
        double rect = max(abs(x),abs(y));
        if( (rect<=rect_max) && (rect>=rect_min) )
        {
          double p = (rect_max-rect)/(rect_max-rect_min);
          double q = 1.-p;
          
          double beta = atan2(y,x);
          double alpha = beta;
          if(y<0)
          {
            alpha = 2.*M_PI+beta;
          }
          
//          cerr << "xy alpha " << x << " " << y << "  :::  " << alpha << " *** " << sin(alpha) << endl;
          double xnew=0.;
          double ynew=0.;
          
//          Choix I
          double r = rect;
          for(int n=0;n<_q[_count].size();n++)
          {
            r += _q[_count][n] * (cos((n+1)*alpha)-1.);
          }
          xnew = q * y  +  p * r*sin(alpha);
          ynew = q * x  +  p * r*cos(alpha);
          nx[i].y() = xnew;
          nx[i].x() = ynew;
        }
      }
      _count++;
    }
};

/*-----------------------------------------*/

class LocalMultiGridMesh : public GascoigneMultiGridMesh
{
  private:
    virtual GascoigneMesh* NewMesh(int dim) 
    {
      if(dim==2)
      {
        return new LocalMesh2d;
      }
      assert(0);
    }

  public:
    LocalMultiGridMesh() {}
    ~LocalMultiGridMesh() {}
};

/*-----------------------------------------*/

class LocalMeshConstructor : public GascoigneMeshConstructor
{
  protected:
    void Construct2d(GascoigneMesh* NM, const LevelMesh2d* LM) const 
    {
      GascoigneMeshConstructor::Construct2d(NM,LM);
      LocalMesh2d* MP = dynamic_cast<LocalMesh2d*>(NM);
      assert(MP);
      MP->SetCoordinates(LM);
    }
    
  public:
    LocalMeshConstructor(const HierarchicalMesh* mm, GascoigneMultiGridMesh* gmg)  : GascoigneMeshConstructor(mm,gmg) {}
    virtual ~LocalMeshConstructor() {}
};

/*-----------------------------------------*/
class LocalMeshAgent : public MeshAgent
{
  protected:
    GascoigneMultiGridMesh* NewMultiGridMesh() { return new LocalMultiGridMesh;}  
    void ReInit() 
    {
      GMG->ReInit(_dimension,HMP->nlevels()-1);
      LocalMeshConstructor MGM(HMP,GMG);
      MGM.BasicInit();
    }
  
  public:
  LocalMeshAgent() : MeshAgent() {}
  ~LocalMeshAgent() {}
  
    void refine_nodes(nvector<int>& refnodes, nvector<int>& coarsenodes)
    {
      cerr << "§§§\n";
      assert(HMP);
      HMP->vertex_patch_refine(refnodes,coarsenodes);
      MeshAgent::ReInit();
    }
    void SetCoordinates() 
    {
      LocalMeshConstructor MGM(HMP,GMG);
      MGM.BasicInit();  
    }
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
