#ifndef  __BenchMarkProblem_h
#define  __BenchMarkProblem_h

#include  "problemdescriptorbase.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include  "navierstokeslps2d.h"

/*-----------------------------------------*/

class BenchMarkDirichletData : public Gascoigne::DirichletData
{
public:

  BenchMarkDirichletData() {}
  std::string GetName() const {return "Bench";}
  void operator()(Gascoigne::DoubleVector& b, const Gascoigne::Vertex2d& v, int color) const 
  {
    double y = v.y();

    b.zero();
    if (color!=80)
      {
        double high = 4.1;
	double vmax = 0.3;
	b[1] = vmax * Gascoigne::ParabelFunction(y,0.,high);
      }
  }
};

/*---------------------------------------------------*/

class BenchMarkProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "NavierStokesBenchmark";}
  void BasicInit(const Gascoigne::ParamFile* pf) 
  {
    GetParamFilePointer()     = pf;
    GetEquationPointer()      = new Gascoigne::NavierStokesLps2d(GetParamFile());
    GetDirichletDataPointer() = new BenchMarkDirichletData;
    
    Gascoigne::ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

#endif
