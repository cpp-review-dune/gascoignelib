
#ifndef  __local_h
#define  __local_h

#include  "ns.h"
#include  "dirichletdata.h"
#include  "problemdescriptorbase.h"
#include  "domainrighthandside.h"
#include  "componentinformationbase.h"
#include <stdio.h>
#include <stdlib.h>
#include "weighteddiracrighthandside.h"
#include "functionals.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

extern double __TIME;
extern double __DT;


class TGV_DD : public DirichletData
{
 protected:
  double __nu;  
 public:
  TGV_DD(const ParamFile* pf)
    {
      DataFormatHandler DFH;
      DFH.insert("nu" ,    &__nu , 0.0);      
      FileScanner FS(DFH, pf, "Equation");
    }
  
  std::string GetName() const {return "TGV-DD";}
  
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();
    double px = v.x()*M_PI;
    double py = v.y()*M_PI;
    
    double F = exp(-2.0*__nu * __TIME * M_PI*M_PI);
    b[1] =  F * (sin(px)*cos(py));
    b[2] = -F * (cos(px)*sin(py));
  }
};


class TaylorGreenVortex : public DomainRightHandSide
{
  double __nu;
  
 public:

  TaylorGreenVortex(const ParamFile* pf)
    {
      DataFormatHandler DFH;
      DFH.insert("nu" ,    &__nu , 0.0);
      FileScanner FS(DFH, pf, "Equation");
      assert(__nu>0);
    }
  

  std::string GetName() const {return "TGV";}
  int GetNcomp() const {return 3;}
    
  virtual double operator()(int c, const Vertex2d& v) const
  {
    double px = v.x()*M_PI;
    double py = v.y()*M_PI;

    double t = 0.0;
    double F = exp(-2.0*__nu * t * M_PI*M_PI);
    if (c==0) return F*F/4.0 * (cos(2.0*px)+cos(2.0*py));
    if (c==1) return  F * (sin(px)*cos(py));
    if (c==2) return -F * (cos(px)*sin(py));
    abort();
  }
};





class MyRHS : public DomainRightHandSide
{
  double __nu;
  
 public:

  MyRHS(const ParamFile* pf)
    {
      DataFormatHandler DFH;
      DFH.insert("nu" ,    &__nu , 0.0);
      FileScanner FS(DFH, pf, "Equation");
      assert(__nu>0);
    }
  

  std::string GetName() const {return "MyDD";}
  int GetNcomp() const {return 3;}
    
  virtual double operator()(int c, const Vertex2d& v) const
  {
    return F(__nu, c,v,__TIME);
  }
};

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
   
    double sc = 1.0;
    if (__TIME<2.0) sc = __TIME/2.0;
    else if (__TIME<6.0) sc = 1.0;
    else sc=(8-__TIME)/2.0;
    sc = __TIME;
    
    if (color==8) 
      b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc; 
  }
};


class MyAdjointDD : public DirichletData
{
 protected:
  
 public:
  MyAdjointDD(const ParamFile* pf)
    {
    }
  
  std::string GetName() const {return "MyAdjointDD";}
  
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();
  }
};

// -----------------------------------------

class ProblemDescriptor2d : public ProblemDescriptorBase
{
 public:
    
  std::string GetName() const {return "primal";}
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new NS(GetParamFile());
    GetRightHandSidePointer() = new MyRHS(GetParamFile());
    GetDirichletDataPointer() = new MyDD(GetParamFile());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

class AdjointProblemDescriptor2d : public ProblemDescriptorBase
{
 public:
    
  std::string GetName() const {return "adjoint";}
  void BasicInit(const ParamFile* pf, DomainRightHandSide& RHS)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new AdjointNS(GetParamFile());
    GetDirichletDataPointer() = new MyAdjointDD(GetParamFile());
    GetRightHandSidePointer() = &RHS;
   
    ProblemDescriptorBase::BasicInit(pf);
  }
  void BasicInit(const ParamFile* pf,  WeightedPointFunctional& WPF)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new AdjointNS(GetParamFile());
    GetDirichletDataPointer() = new MyAdjointDD(GetParamFile());
    GetRightHandSidePointer() = new WeightedDiracRightHandSide();
    dynamic_cast< WeightedDiracRightHandSide*>(GetRightHandSidePointer())->BasicInit(&WPF);
    
   
    ProblemDescriptorBase::BasicInit(pf);
  }
  void BasicInit(const ParamFile* pf, BoundaryRightHandSide& WPF)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new AdjointNS(GetParamFile());
    GetDirichletDataPointer() = new MyAdjointDD(GetParamFile());
    GetBoundaryRightHandSidePointer() = &WPF;

    
   
    ProblemDescriptorBase::BasicInit(pf);
  }
};




#endif
