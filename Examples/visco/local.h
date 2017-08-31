
#ifndef  __local_h
#define  __local_h

#include  "fsi.h"
#include  "dirichletdata.h"
#include  "problemdescriptorbase.h"
#include  "domainrighthandside.h"
#include  "componentinformationbase.h"
#include <stdio.h>
#include <stdlib.h>


using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

extern double __TIME;



class MyDD : public DirichletData
{
 protected:
  double vmean;
  
 public:
  MyDD(const ParamFile* pf)
    {
      DataFormatHandler DFH;
      DFH.insert("vmean" ,    &vmean , 0.0);
      FileScanner FS(DFH, pf, "Equation");
    }
  
  std::string GetName() const {return "MyDD";}
  
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();
    if (color==8)
      b[1] = 8.0 * (1.0+ tanh( 8.0 * (__TIME - 0.5))) * v.x()*v.x()*(1.0-v.x())*(1.0-v.x());
  }
};
class MyStress : public DirichletData
{
 protected:
  
 public:

  
  std::string GetName() const {return "MyStress";}
  
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b[0]=1.0;
    b[2]=1.0;
  }
};


// -----------------------------------------

class VelProblem : public ProblemDescriptorBase
{
 public:
    
  std::string GetName() const {return "vel";}
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer()      = new VelEQ(GetParamFile());
    GetDirichletDataPointer() = new MyDD(GetParamFile());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

class StressProblem : public ProblemDescriptorBase
{
 public:
    
  std::string GetName() const {return "stress";}
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer()      = new StressEQ(GetParamFile());
    GetDirichletDataPointer() = new MyStress;
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};





#endif
