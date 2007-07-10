#ifndef  __local_h
#define  __local_h

#include  "stdtimeloop.h"
#include  "navierstokesgls2d.h"
#include  "dirichletdata.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include  "problemdescriptorbase.h"
#include  "navierstokeslps2d.h"
#include  "stokeslps2d.h"
#include  "hierarchicalmesh2d.h"
#include  "nsface.h"

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
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const {

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
  void BasicInit(const ParamFile* pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer()      = new NavierStokes2d(GetParamFile());
    GetFaceEquationPointer()  = new NavierStokesFace2d(GetParamFile());
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

      AddShape(80,&RK);
    }
  void BasicInit(const ParamFile* paramfile)
    {
      std::string inpname("nsbench4.inp");
      
      HMP = new HierarchicalMesh2d;
      map<int,BoundaryFunction<2>* >::const_iterator p;
      for(p=GetShapes2d().begin();p!=GetShapes2d().end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
      int patchdepth = 1;
      int etapatcher = 1;
      HMP->SetParameters(inpname,patchdepth,etapatcher);
      int prerefine  = 1;
      HMP->global_refine(prerefine);
      
      GMG = NewMultiGridMesh();
      
      ReInit();
    }
};

/* ----------------------------------------- */

#include  "domainfunctional.h"
#include  "residualfunctional.h" 
#include  "dirichletdatabycolor.h"

/* class LocalDomainFunctionals_FlagForce : public virtual Gascoigne::ResidualFunctional */
/* { */
/*  public: */
/*   std::string _s_force_type; */
/*   LocalDomainFunctionals_FlagForce(std::string s_force_type) : ResidualFunctional() { */
/*     _s_force_type = s_force_type; */

/*     if(s_force_type == "drag") __comps.push_back(1); */
/*     if(s_force_type == "lift") __comps.push_back(2); */
/*     assert(__comps.size()==1); */

/*     __cols.insert(80); */
/*     __scales.push_back(1); */
/*     ExactValue() = 0.; */

/*     __DD  = new Gascoigne::DirichletDataByColor(GetComps(),GetColors(),GetScales()); */
/*   } */

/*   std::string GetName() const { */
/*     return "LocalDomainFunctionals_FlagForce:"+_s_force_type; */
/*   } */
/* }; */

/* ----------------------------------------- */


#endif