#ifndef  __local_h
#define  __local_h

#include  "stdloop.h"
#include  "navierstokesgls3d.h"
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
    vmax = 0.45;
  }
  std::string GetName() const {return "Bench";}
  void operator()(DoubleVector& b, const Vertex3d& v, int color) const {

    double x = v.x();  double y = v.y(); double z = v.z();
    
    b.zero();
    if (x==0.)
      {
	double high = 4.1;
	b[1] = vmax * ParabelFunction(y,0.,high) * ParabelFunction(z,0.,high);
      }
  }
};

/* ----------------------------------------- */

class ProblemDescriptor : public ProblemDescriptorBase
{
public:
    
    std::string GetName() const {return "Local";}
    void BasicInit(const ParamFile* pf) {
      GetParamFilePointer() = pf;
      GetEquationPointer() = new NavierStokesGls3d(GetParamFile());
      GetDirichletDataPointer() = new BenchMarkDirichletData();
      ProblemDescriptorBase::BasicInit(pf);
    }
};


/* ----------------------------------------- */

class RunderKreis : public BoundaryFunction<2>
{
  double squareradius;
  Vertex2d center;
public :

  std::string GetName() const { return "RunderKreis";}

  void BasicInit(Vertex2d c, double r) 
    {
      center = c; 
      squareradius = r;
    }
  double operator()(const Vertex2d& c) const 
    {
      double r = - squareradius;
      for (int i=0; i<2; i++)
	{
	  double dx = c[i]-center[i];
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
      AddShape(80,&RK);
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
