/*----------------------------   equation.h     ---------------------------*/
/*      $Id: myequation.h,v 1.1 2009/10/26 08:42:26 richter Exp $                 */
#ifndef __equation_H
#define __equation_H
/*----------------------------   equation.h     ---------------------------*/


#define  __MyEquation_h

#include  "paramfile.h"
#include  "equation.h"

using namespace std;

/*-----------------------------------------*/
 extern bool FIRSTDUAL;
extern bool LASTDUAL;
namespace Gascoigne
{
  
  class MyEquation : public virtual Equation
  {
  protected:

    mutable FemFunction* oldu;
    double epsilon;
  public:
    
    MyEquation() { abort(); }
    MyEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldu") != q.end() );
      oldu = &q["oldu"];
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  };


  class MyDualEquation : public virtual Equation
  {
  protected:

     mutable FemFunction* oldz;
     mutable FemFunction* u1;
     mutable FemFunction* u2;
     mutable FemFunction* u3;
   
     double epsilon;
  public:
    
    MyDualEquation() { abort(); }
    MyDualEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyDualEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldz") != q.end() ); 
      oldz = &q["oldz"];
      
      if (!LASTDUAL)
	{
	  assert(q.find("u1") != q.end() ); 
	  u1 = &q["u1"];
	}
      else
	u1 = NULL;
      
      assert(q.find("u2") != q.end() ); 
      u2 = &q["u2"];

      if (!FIRSTDUAL)
	{
	  assert(q.find("u3") != q.end() ); 
	  u3 = &q["u3"];
	}
      else
	u3 = NULL;
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w) const;
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Nonlinear_Matrix(EntryMatrix&A, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction &M,const TestFunction& N,double w) const;
    void Matrix(EntryMatrix& A, const FemFunction& Z, const TestFunction& M, const TestFunction& N) const;

  };


 
}


/*----------------------------   equation.h     ---------------------------*/
/* end of #ifndef __equation_H */
#endif
/*----------------------------   equation.h     ---------------------------*/
