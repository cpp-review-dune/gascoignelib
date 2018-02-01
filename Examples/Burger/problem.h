/*----------------------------   problem.h     ---------------------------*/
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
  
  class MyRhs : public virtual DomainRightHandSide
  {
    
  public:
    
  MyRhs(const ParamFile* pf) : DomainRightHandSide()
      {
      }
    
    std::string GetName() const {return "MyRhs";}    
    int GetNcomp() const {return 2; }
    
    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    {
      double x=( v.x()+0.5)*( v.x()+0.5)+(v.y()-0.2)*(v.y()-0.2);
      double y= ( v.x()-0.5)*( v.x()-0.5)+(v.y()+0.2)*(v.y()+0.2);
    
      double w=0.0;
      w=sin(TIME);
      
      
      b[0] +=w*(exp(-2*x)-exp(-2*y)) * N.m(); 
      b[1] +=w*(exp(-2*y)-exp(-2*x)) * N.m();

      
	
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
     	
    }
    
  };
  
  class MyDualRhsDom : public DomainRightHandSide
  {
    
    
  public:
    
    int GetNcomp()        const { return 2; }
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
	   //if(v.x()>0.0 && v.x()<0.5 &&v.y()>0.0 && v.y()<0.5) 
	     //{
             int x=v.x()-0.5;
             int y =v.y()-0.5;
	       b[0]+= N.x()*exp(-x*x-y*y);
	       b[1]+=-N.y()*exp(-x*x-y*y);
	  //   }
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

  class DualProblem : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "Dual  Problem";}
    void BasicInit(const ParamFile* pf)
    {
      // equation to solve
      GetEquationPointer()      = new MyDualEquation(pf); 
      GetRightHandSidePointer() = new MyDualRhsDom(pf);
      GetDirichletDataPointer() = new MyDualDD();
      // 
      ProblemDescriptorBase::BasicInit(pf);
    }
  };


  
  
}



/*----------------------------   problem.h     ---------------------------*/
/* end of #ifndef __problem_H */
#endif
/*----------------------------   problem.h     ---------------------------*/
  
