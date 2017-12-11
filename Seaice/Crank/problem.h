#/*----------------------------   problem.h     ---------------------------*/
/*      $Id: problem.h,v 1.2 2009/10/26 09:07:16 richter Exp $                 */
#ifndef __problem_H
#define __problem_H
/*----------------------------   problem.h     ---------------------------*/


#include  "problemdescriptorbase.h"
#include  "myequation.h"
#include  "dirichletdata.h"
#include  "domainrighthandside.h"
#include  "zerodirichletdata.h"



extern double TIME,DT,DTM;
using namespace std;
namespace Gascoigne
{
  


      
  
  class MyNonLinearRHS : public virtual DomainRightHandSide
  {

    mutable FemFunction *Pu_kM;

    void SetFemData(FemData& q) const 
    { 
      assert(q.find("Pu_kM")!=q.end()); 
      Pu_kM = &q["Pu_kM"];

    }
    
  public:
    
  MyNonLinearRHS(const ParamFile* pf) : DomainRightHandSide()
      {
	DataFormatHandler DFH;
	FileScanner FS(DFH,pf,"Equation");
	FS.NoComplain();
      }
    
    std::string GetName() const {return "MyRhs";}    
    int GetNcomp() const {return 1; }
    
    void operator()(VectorIterator b, const TestFunction& N,const TestFunction &M, const Vertex2d& v) const 
    {
      b[0] += ((*Pu_kM)[0].x()*(*Pu_kM)[0].x()+(*Pu_kM)[0].y()*(*Pu_kM)[0].y() )* N.m() ;
    }
  };

  

  class MyRhs : public virtual DomainRightHandSide
  {
    
  public:
    
  MyRhs(const ParamFile* pf) : DomainRightHandSide()
      {
      }
    
    std::string GetName() const {return "MyRhs";}    
    int GetNcomp() const {return 1; }
    
    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    {
      double x= v.x()-0.75;
      double y= v.y()-0.5;
      b[0] += N.m() * exp(-10.0*(x*x+y*y));
      return;
      

      
      /*
      double a=50;
      double  cg0 = -pow(0.1e1 + 0.50e2 * pow(x - 0.1e1 / 0.2e1 - cos(0.2e1 * 0.3141592654e1 * t) / 0.4e1, 0.2e1) + 0.50e2 * pow(y - 0.1e1 / 0.2e1 - sin(0.2e1 * 0.3141592654e1 * t) / 0.4e1, 0.2e1), -0.2e1) * (0.50e2 * (x - 0.1e1 / 0.2e1 - cos(0.2e1 * 0.3141592654e1 * t) / 0.4e1) * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * t) - 0.50e2 * (y - 0.1e1 / 0.2e1 - sin(0.2e1 * 0.3141592654e1 * t) / 0.4e1) * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * t)) - 0.2e1 * pow(0.1e1 + 0.50e2 * pow(x - 0.1e1 / 0.2e1 - cos(0.2e1 * 0.3141592654e1 * t) / 0.4e1, 0.2e1) + 0.50e2 * pow(y - 0.1e1 / 0.2e1 - sin(0.2e1 * 0.3141592654e1 * t) / 0.4e1, 0.2e1), -0.3e1) * pow(0.100e3 * x - 0.50e2 - 0.25e2 * cos(0.2e1 * 0.3141592654e1 * t), 0.2e1) + 0.200e3 * pow(0.1e1 + 0.50e2 * pow(x - 0.1e1 / 0.2e1 - cos(0.2e1 * 0.3141592654e1 * t) / 0.4e1, 0.2e1) + 0.50e2 * pow(y - 0.1e1 / 0.2e1 - sin(0.2e1 * 0.3141592654e1 * t) / 0.4e1, 0.2e1), -0.2e1) - 0.2e1 * pow(0.1e1 + 0.50e2 * pow(x - 0.1e1 / 0.2e1 - cos(0.2e1 * 0.3141592654e1 * t) / 0.4e1, 0.2e1) + 0.50e2 * pow(y - 0.1e1 / 0.2e1 - sin(0.2e1 * 0.3141592654e1 * t) / 0.4e1, 0.2e1), -0.3e1) * pow(0.100e3 * y - 0.50e2 - 0.25e2 * sin(0.2e1 * 0.3141592654e1 * t), 0.2e1);  
      double t2 = 0.2e1 * 0.3141592654e1 * t;
      double t3 = cos(t2);
      double t5 = x - 0.1e1 / 0.2e1 - t3 / 0.4e1;
      double t6 = t5 * t5;
      double t8 = sin(t2);
      double t10 = y - 0.1e1 / 0.2e1 - t8 / 0.4e1;
      double t11 = t10 * t10;

      double t13 = 0.1e1 + 0.50e2 * t6 + 0.50e2 * t11;
      double t14 = t13 * t13;
      double t15 = 0.1e1 / t14;
      double t23 = 0.1e1 / t14 / t13;
      double t27 = pow(0.100e3 * x - 0.50e2 - 0.25e2 * t3, 0.2e1);
      double t34 = pow(0.100e3 * y - 0.50e2 - 0.25e2 * t8, 0.2e1);
      double t37 = -0.50e2 * (-t3 * 0.3141592654e1 * t10 + t8 * 0.3141592654e1 * t5) * t15 - 0.2e1 * t27 * t23 + 0.200e3 * t15 - 0.2e1 * t34 * t23;
 
      b[0] +=t37*N.m();
      */
    }
  };

  
   
       
   
  class MyDD : virtual public DirichletData
  {
  public:
    std::string GetName() const {return "DD";}  
  
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const 
    {
      b.zero();
      return;
      /*
      double x=v.x();
      double y=v.y();
      double t=0.0*TIME;
      double a=50;
      double 
	cg1 = 0.1e1 / (0.1e1 + 0.50e2 * pow(x - 0.1e1 / 0.2e1 - cos(0.2e1 * M_PI* t) / 0.4e1, 0.2e1) + 0.50e2 * pow(y - 0.1e1 / 0.2e1 - sin(0.2e1 * M_PI * t) / 0.4e1, 0.2e1)); 
      b[0]=cg1;
      // b[0]=0.0;	
      */
    }
    
  };
  
  class MyDualDD : virtual public DirichletData
  {
  public:
    std::string GetName() const {return "DD";}  
  
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const 
    {
      b.zero();
    }
    
  };

  
  class MyDualRhs : public BoundaryRightHandSide
  {
    
  public:
    int GetNcomp()        const { return 1; }
    std::string GetName() const { return "Dual Rhs"; }
        
  MyDualRhs(const ParamFile* paramfile) : BoundaryRightHandSide()
      {
	DataFormatHandler DFH;
	FileScanner FS(DFH);
	FS.NoComplain();
	FS.readfile(paramfile, "Equation");
      }


    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v, const Vertex2d& n, int color) const
    {
      double c= 0.249702115724673;
      
     if(v.y()==0)
	{
	  b[0]+= DTM*1/c * N.y();
	}
    }
  };
  
class MyDualRhsDom : public DomainRightHandSide
  {
    
  public:
    

   

    int GetNcomp()        const { return 1; }
    std::string GetName() const { return "Dual Rhs"; }
    
    
    
  MyDualRhsDom(const ParamFile* paramfile) : DomainRightHandSide()
      {
	DataFormatHandler DFH;
	FileScanner FS(DFH);
	FS.NoComplain();
	FS.readfile(paramfile, "Equation");
      }

    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
        {
	  
	  
	  double x=v.x()-0.25;
          double  y=v.y()-0.25;
	  
	  b[0]+=DTM*exp(-10*(x*x+y*y))*N.m();
	    
	}
  };
 





  // main class for defining the problem to solve
  class LaplaceTProblem : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "LaplaceTProblem";}
    void BasicInit(const ParamFile* pf)
    {
      // equation to solve
      GetEquationPointer()      = new MyEquation(pf);

      GetRightHandSidePointer() = new MyRhs(pf);

      GetDirichletDataPointer() = new MyDD();
      // 
      ProblemDescriptorBase::BasicInit(pf);
    }
  };
  


  // main class for defining the problem to solve
  class DualProblem : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "Dual  Problem";}
    void BasicInit(const ParamFile* pf)
    {
      // equation to solve
      GetEquationPointer()      = new MyDualEquation(pf); // spaeter dualequation, ohne spruenge ddd=1

      // definition for right hand side
      // GetBoundaryRightHandSidePointer() = new MyDualRhs(pf);
      GetRightHandSidePointer() = new MyDualRhsDom(pf);
      GetDirichletDataPointer() = new MyDualDD();
      // 
      ProblemDescriptorBase::BasicInit(pf);
    }
  };

   // main class for defining the problem to solve
  class DualHelpProblem : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "Dual Help  Problem";}
    void BasicInit(const ParamFile* pf)
    {
      // equation to solve
      GetEquationPointer()      = new MyDualHelpEquation(pf); 

      // definition for right hand side
      GetRightHandSidePointer() = new  MyRhs(pf);

      GetDirichletDataPointer() = new MyDD();
      // 
      ProblemDescriptorBase::BasicInit(pf);
    }
  };
  
 /* class NonlinearProblem : public ProblemDescriptorBase */
 /*  { */
 /*  public: */
  
 /*    std::string GetName() const {return "Nonlinear Problem";} */
 /*    void BasicInit(const ParamFile* pf) */
 /*    { */
 /*      // equation to solve */
 /*      GetEquationPointer()      = new MyNonLinearEquation(pf);  */

 /*      // definition for right hand side */
 /*      GetRightHandSidePointer() = new  MyNonLinearRHS(pf); */

 /*      GetDirichletDataPointer() = new MyDD(); */
 /*      //  */
 /*      ProblemDescriptorBase::BasicInit(pf); */
 /*    } */
 /*  }; */
  


}



/*----------------------------   problem.h     ---------------------------*/
/* end of #ifndef __problem_H */
#endif
/*----------------------------   problem.h     ---------------------------*/
  
