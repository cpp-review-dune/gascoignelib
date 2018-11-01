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
  //    double x=v.x()-0.75;
  //    double y=v.y()-0.25;
    
      
      // b[0] +=(exp(-10*x*x-10*y*y)) * N.m(); 
       // b[1] +=0;
    }
  };

  
   
       
   
  class MyDD : virtual public DirichletData
  {
  public:
    std::string GetName() const {return "DD";}  
  
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const 
    {
        double x=v.x();
        double y=v.y();
        double x1=x-0.5;
        
      b[0]=(exp(-10*x1*x1-10*y*y))*(1-x*x)*(1-y*y); 
      b[1]=(exp(-10*x1*x1-10*y*y))*(1-x*x)*(1-y*y); 
     	
    }
    
  };
  
  class MyDualRhsDom : public DomainRightHandSide
  {
    
    mutable FemFunction *H;
    
      
      
  public:
      
    void SetFemData(FemData& q) const 
    {
      assert(q.find("H") != q.end() ); 
      H= &q["H"];
    }
    
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
      double x = v.x()-0.75;
      double y = v.y()-0.25;
      
      b[0]+= N.m()*(exp(-10*x*x-10*y*y))*(1-x*x)*(1-y*y);
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
  
  
  // Burger gleichung
   
  class MyBurgerRhs : public virtual DomainRightHandSide
  {
     
  public:
     
  MyBurgerRhs(const ParamFile* pf) : DomainRightHandSide()
      {
      }
     
    std::string GetName() const {return "MyRhs";}    
    int GetNcomp() const {return 2; }
    
    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    {
        double x = v.x();
      double y = v.y(); 
      
     b[0]+= -v.y()*N.m()*(1-x*x)*(1-y*y) ; 
     b[1]+=  v.x()*N.m()*(1-x*x)*(1-y*y) ;
     
      
    }
  };

  
   
       
   
  class MyBurgerDD : virtual public DirichletData
  {
  public:
    std::string GetName() const {return "DD";}  
  
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const 
    {
      b.zero();
      return;
     	
    }
    
  };
  
  class MyBurgerDualRhsDom : public DomainRightHandSide
  {
    
    mutable FemFunction *V;
  public:
    
    int GetNcomp()        const { return 2; }
    std::string GetName() const { return "Dual Rhs"; }
    
    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("V") != q.end() ); 
      V = &q["V"];
      }

    
  MyBurgerDualRhsDom(const ParamFile* paramfile) : DomainRightHandSide()
      {
	DataFormatHandler DFH;
	FileScanner FS(DFH);
	FS.NoComplain();
	FS.readfile(paramfile, "Equation");
      }

    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    {
        
        /*
      b[0] +=  ((*V)[0].x() - (*V)[1].y()) * N.x();
      b[1] += -((*V)[0].x() - (*V)[1].y()) * N.y();
      
      b[0] += ((*V)[1].y()+(*V)[0].x()) * N.x();
      b[1] += ((*V)[1].y()+(*V)[0].x()) * N.y();
      */
      /* return 0.5 * (pow(U[0].x()-U[1].y(),2.0) + */
      /* 		    pow(U[1].y()+U[0].x(),2.0)); */
      
    }
  };
  


  class MyBurgerDualDD : virtual public DirichletData
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
      GetEquationPointer()      = new TransportEquation(pf);
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
  
  
  // Burger---
  
  class BurgerProblem : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "BurgerProblem";}
    void BasicInit(const ParamFile* pf)
    {
      // equation to solve
      GetEquationPointer()      = new MyBurgerEquation(pf);
      GetRightHandSidePointer() = new MyBurgerRhs(pf);
      GetDirichletDataPointer() = new MyBurgerDD();
      // 
      ProblemDescriptorBase::BasicInit(pf);
    }
  };

  class BurgerDualProblem : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "BurgerDualProblem";}
    void BasicInit(const ParamFile* pf)
    {
      // equation to solve
      GetEquationPointer()      = new MyBurgerDualEquation(pf); 
     GetRightHandSidePointer() = new MyBurgerDualRhsDom(pf);
      GetDirichletDataPointer() = new MyBurgerDualDD();
      // 
      ProblemDescriptorBase::BasicInit(pf);
    }
  };
  
  
  
}



/*----------------------------   problem.h     ---------------------------*/
/* end of #ifndef __problem_H */
#endif
/*----------------------------   problem.h     ---------------------------*/
  
