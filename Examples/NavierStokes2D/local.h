#ifndef  __local_h
#define  __local_h

#include  "stdloop.h"
#include  "navierstokesgls2d.h"
#include  "dirichletdata.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"

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

    double x = v.x();  double y = v.y();
    
    b.zero();
    if (color!=80)
      {
	double high = 4.1;
	b[1] = vmax * ParabelFunction(y,0.,high);
      }
    if(b.size()>3)
      {
	double tmin = 300.;
	double tmax = 600.;
	//double tmax = 300.;
	
	double x0 = 2.;
	double y0 = 2.;
	
	if ((x-x0)*(x-x0)+(y-y0)*(y-y0) > 1. ) b[3]= tmin;
	else b[3] = tmax;
      }
    if(b.size()>4)
      {
	double T = b[3];
	b[4] = 273./ T;
      }
  }
};

/* ----------------------------------------- */

class ProblemDescriptor1 : public ProblemDescriptorInterface
{
protected:

  void ConstructEquation() {
    GetEquationPointer() = new NavierStokesGls2d(GetParamFile());
  }
  void ConstructDirichletData() {
    GetDirichletDataPointer() = new BenchMarkDirichletData();
  }
  
public:
    
    std::string GetName() const {return "Local";}
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
