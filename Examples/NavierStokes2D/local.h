#ifndef  __local_h
#define  __local_h

#include  "stdloop.h"
#include  "navierstokesgls2d.h"
#include  "dirichletdata.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include  "problemdescriptorbase.h"

using namespace std;
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
  void operator()(Vector& b, const Vertex2d& v, int color) const {

    double y = v.y();
    
    b.zero();
    if (color!=80)
      {
	double high = 4.1;
	b[1] = vmax * ParabelFunction(y,0.,high);
      }
  }
};

/* ----------------------------------------- */

class ProblemDescriptor : public ProblemDescriptorBase
{
public:
    
    std::string GetName() const {return "Local";}
    void BasicInit(const Gascoigne::ParamFile* pf) {
      GetEquationPointer() = new NavierStokesGls2d(GetParamFile());
      GetDirichletDataPointer() = new BenchMarkDirichletData();
      ProblemDescriptorBase::BasicInit(pf);
    }
};


/* ----------------------------------------- */

class RunderKreis : public BoundaryFunction<2>
{
  double   _r;
  Vertex2d _c;
public :

  std::string GetName() const { return "RunderKreis";}
  void BasicInit(Vertex2d c, double r) {
      _c = c; 
      _r = r;
    }
  double operator()(const Vertex2d& c) const {
      double r = - _r;
      for (int i=0; i<2; i++)
	{
	  double dx = c[i]-_c[i];
	  r += dx * dx;
	}
      return r;
    }
};

/*---------------------------------------------------*/

class BenchMarkMeshAgent : public MeshAgent
{
 protected:
  
  RunderKreis RK;

 public:
  
  BenchMarkMeshAgent() : MeshAgent()
    {
      double r = 0.25;
      Vertex2d v(2.,2.);
      RK.BasicInit(v,r);

      int dim=2;
      int prerefine=1;
      std::string inpname("nsbench4.inp");

      AddShape(80,&RK);
      BasicInit(dim,inpname,prerefine);
    }
};

/* ----------------------------------------- */

class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile) 
    {
      GetMeshAgentPointer() = new BenchMarkMeshAgent;
      StdLoop::BasicInit(paramfile);
    }
};


#endif
