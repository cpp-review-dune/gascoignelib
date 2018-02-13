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
    mutable FemFunction* H;
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
       assert(q.find("H") != q.end() );
      H = &q["H"];
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
     
    mutable FemFunction* h;
    mutable FemFunction* w;
   
    mutable FemFunction* newH;
    mutable FemFunction* oldW;
    
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
      
      assert(q.find("w") != q.end() ); 
	  w = &q["w"];
      
	}
      else{
	u1 = NULL;
    w=NULL;
          
    
    }
    
      assert(q.find("u2") != q.end() ); 
      u2 = &q["u2"];
        assert(q.find("h") != q.end() ); 
	  h = &q["h"];
      assert(q.find("newH") != q.end() ); 
	  newH = &q["newH"];

      if (!FIRSTDUAL)
	{
	  assert(q.find("u3") != q.end() ); 
	  u3 = &q["u3"];
      
      assert(q.find("oldW") != q.end() ); 
	  oldW = &q["oldW"];
      
      
	}
      else{
	u3 = NULL;
    oldW=NULL;
    newH=NULL;
      }
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const;
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Nonlinear_Matrix(EntryMatrix&A, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction &M,const TestFunction& N,double w,int DTM) const;
    void Matrix(EntryMatrix& A, const FemFunction& Z, const TestFunction& M, const TestFunction& N) const;

  };



    class MyTransportEquation : public virtual Equation
  {
  protected:

    mutable FemFunction* oldh;
    mutable FemFunction* V;
    // mutable FemFunction V;
    double epsilon;
  public:
    
    MyTransportEquation() { abort(); }
    MyTransportEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyTransortEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldh") != q.end() );
      oldh = &q["oldh"];
      assert(q.find("V") != q.end() );
      V = &q["V"];
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  };
 
}


/*----------------------------   equation.h     ---------------------------*/
/* end of #ifndef __equation_H */
#endif
/*----------------------------   equation.h     ---------------------------*/
