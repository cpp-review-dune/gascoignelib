#ifndef  __local_h
#define  __local_h

#include  "problemdescriptorbase.h"
#include  "localequation.h"
#include  "dirichletdata.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include  "stdloop.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

class LocalDirichletData : public DirichletData
{
 public:
  std::string GetName() const {return "Local";}
  void operator()(DoubleVector& b, const Vertex2d& v, int col) const {
    if(col!=80) {
      b[0] = 0.;
    } else {
      b[0] = 1.;
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
      GetEquationPointer() = new LocalEquation(GetParamFile());
      GetDirichletDataPointer() = new LocalDirichletData;
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
  void run(const ProblemDescriptorInterface* PD);
};


#endif
