/*----------------------------   problem.h     ---------------------------*/
/*      $Id: problem.h,v 1.4 2013/05/08 12:51:16 richter Exp $                 */
#ifndef __problem_H
#define __problem_H
/*----------------------------   problem.h     ---------------------------*/

#include  "problemdescriptorbase.h"
#include  "myequation.h"
#include  "dirichletdata.h"
#include  "domainrighthandside.h"
#include  "domainmeanfunctional.h"
#include  "weighteddiracrighthandside.h"

/*---------------------------------------------------*/

extern double __GLOBAL_TIME;
extern double __GLOBAL_DT;

namespace Gascoigne
{
  
// --------------------------------------------------
  
  static double scale_t(double x)
  {
    return 1.0;
  }
  
  class MyFunc : public DomainFunctional
  {
  public:
    std::string GetName() const {return "MyExact";}
    double J(const FemFunction& U, const Vertex2d& v) const
    {
      double sc = exp(- M_PI*M_PI * (v.x()*v.x() + (v.y()+0.5)*(v.y()+0.5)));
      return sc * U[0].m();
      
      /* double t= __GLOBAL_TIME; */
      /* double X0 = 0.0; */
      /* double Y0 = -0.5; */
      /* double sc = 0.5 + 0.5*erf(5.0*(t-M_PI)); */
      /* sc = 1.0; */
      /* if (t<M_PI*0.5) sc = 0.0; */
      
      
      /* return sc * (U[0].x()*U[0].x() + U[0].y()*U[0].y()); */
    }
  };

  /* class MyFuncRhs : public WeightedDiracRightHandSide */
  /*   { */
  /*     WeightedPointFunctional F; */
  /*   public: */
  /*     MyFuncRhs() */
  /* 	{ */
  /* 	  std::vector<Vertex2d> vv;vv.push_back(Vertex2d(0.0,-0.5)); */
  /* 	  std::vector<double>   ww;ww.push_back(1.0); */
  /* 	  std::vector<int>      cc;cc.push_back(0); */
  /* 	  F.BasicInit(vv,cc,ww); */
  /* 	  BasicInit(&F); */
  /* 	} */
 
  /*     std::string GetName() const {return "MyFuncRHS";} */
  /*     int GetNcomp() const { return 1; } */
  /*   }; */
  class MyFuncRhs : public DomainRightHandSide
    {
      mutable FemFunction *__U;
    public:
      MyFuncRhs() { __U=0; }
      void SetFemData(FemData& q) const
      { assert(q.find("u")!=q.end()); __U = &q["u"]; }
 
      std::string GetName() const {return "MyFuncRHS";}
      int GetNcomp() const { return 1; }
      void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const
      {
	double sc = exp(- M_PI*M_PI * (v.x()*v.x() + (v.y()+0.5)*(v.y()+0.5)));
	
	b[0] += sc * N.m();
	
	/* double t= __GLOBAL_TIME; */
	/* double X0 = 0.0; */
	/* double Y0 = -0.5; */

	/* double sc = 0.5 + 0.5*erf(5.0*(t-M_PI)); */

	/* sc = 1.0; */
	/* if (t<M_PI*0.5) sc = 0.0; */

	


	/* b[0] += sc * (2.0 * (*__U)[0].x() * N.x() + */
	/* 	      2.0 * (*__U)[0].y() * N.y()); */
      }
    };
  
  class MyInitial : public DomainRightHandSide
  {
  public:
    
    std::string GetName() const {return "MyInitial";}
    int GetNcomp() const { return 1; }
    double operator()(int c, const Gascoigne::Vertex2d& v) const
    {      
      /* return (v.x()*v.x()-1.0)*(v.y()*v.y()-1.0) * exp(-(v.x()-0.5)*(v.x()-0.5)*5.0 - v.y()*v.y()*5.0); */


      /* if ((fabs(v.x()-0.0)<0.25)&& */
      /* 	  (fabs(v.y()-0.5)<0.25)) return 1.0; */
      /* return 0.0; */

      
      return (v.x()*v.x()-1)*(v.y()*v.y()-1) *
      	exp(-M_PI*M_PI * (pow(v.x(),2.0) + pow (v.y()-0.5,2.0)));
    }
  };
  
  class MyInitial2 : public WeightedDiracRightHandSide
  {
    WeightedPointFunctional F;
  public:
    MyInitial2()
      {
	std::vector<Vertex2d> vv;vv.push_back(Vertex2d(0.0,0.5));
	std::vector<double>   ww;ww.push_back(1.0);
	std::vector<int>      cc;cc.push_back(0);
	F.BasicInit(vv,cc,ww);
	BasicInit(&F);
      }
 
    std::string GetName() const {return "MyInitial2";}
    int GetNcomp() const { return 1; }
  };


  class MyRHS : public Gascoigne::DomainRightHandSide
    {
    public:
      std::string GetName() const {return "MyDD";}
      int GetNcomp() const { return 1; }
      double operator()(int c, const Gascoigne::Vertex2d& v) const
      {
	return 0.0;
      }
    };
  
  
  class MyDD : public Gascoigne::DirichletData
    {
    public:

    public:
      std::string GetName() const {return "MyDD";}
      
      void operator()(Gascoigne::DoubleVector& b, const Gascoigne::Vertex2d& v, int col) const
      {
	b.zero();
      }
    };
  
  /*---------------------------------------------------*/
  
  class ProblemDescriptor : public Gascoigne::ProblemDescriptorBase
    {
    public:
      std::string GetName() const {return "Test";}
      void BasicInit(const Gascoigne::ParamFile* pf)
      {
	GetEquationPointer()      = new MyEquation(pf);
	GetRightHandSidePointer() = new MyRHS;
	GetDirichletDataPointer() = new MyDD;
	GetInitialConditionPointer() = new MyInitial;
	ProblemDescriptorBase::BasicInit(pf);
      }
    };
}



/*----------------------------   problem.h     ---------------------------*/
/* end of #ifndef __problem_H */
#endif
/*----------------------------   problem.h     ---------------------------*/
