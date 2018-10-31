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

      
    mutable FemFunction* u0;  
    mutable FemFunction* oldu;
    mutable FemFunction* h;
     mutable FemFunction* oldh;
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
      
      assert(q.find("u0") != q.end() );
      u0 = &q["u0"];
       assert(q.find("h") != q.end() );
      h = &q["h"];
      assert(q.find("oldh") != q.end() );
      oldh = &q["oldh"];
      
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
    void Zeit(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const;
   void Zeit_Matrix(EntryMatrix&A, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z,const TestFunction &M, const TestFunction& N,double w,int DTM) const;

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
    
    mutable FemFunction* oldh;
    mutable FemFunction* oldoldz;
    mutable FemFunction* Hnewnew;
    
    
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
     assert(q.find("u2") != q.end() ); 
      u2 = &q["u2"];
        assert(q.find("h") != q.end() ); 
	  h = &q["h"]; 
       assert(q.find("oldh") != q.end() ); 
	  oldh = &q["oldh"];
      
      
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
    
      
     

      if (!FIRSTDUAL)
	{
	  assert(q.find("u3") != q.end() ); 
	  u3 = &q["u3"];
      
      assert(q.find("oldW") != q.end() ); 
	  oldW = &q["oldW"];
      
      
       assert(q.find("newH") != q.end() ); 
	  newH = &q["newH"];
      
      
	}
      else{
	u3 = NULL;
    oldW=NULL;
    newH=NULL;
    oldoldz=NULL;
    Hnewnew=NULL;
      }
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const;
    
    void Kopplung(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const;
    
  
    
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

     class MyDualTransportEquation : public virtual Equation
  {
  protected:

    mutable FemFunction* oldw;
   
 
    mutable FemFunction* u1;
    mutable FemFunction* u2;
    mutable FemFunction* u3;
    mutable FemFunction* oldz;
    mutable FemFunction* z;
    
    double epsilon;
  public:
    
    MyDualTransportEquation() { abort(); }
    MyDualTransportEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyDualTransortEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldw") != q.end() );
      oldw = &q["oldw"];
      assert(q.find("u2") != q.end() );
      u2 = &q["u2"];
      if (!LASTDUAL){
        assert(q.find("u1") != q.end() );
       u1 = &q["u1"];  
       assert(q.find("oldz") != q.end() );
      oldz= &q["oldz"];    
      
    }
      else
      {
      u1=NULL;    
	
    oldz=NULL;
      } 
      
      if(!FIRSTDUAL){
          
    assert(q.find("u3") != q.end() );
      u3 = &q["u3"];
       assert(q.find("z") != q.end() );
      z= &q["z"];   
      }
      else{
          
    u3=NULL;
    z=NULL;   
    }
     
      
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
   
    void Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const;
    void Nonlinear_Matrix(EntryMatrix&A, double s, const FemFunction &U1, const FemFunction& U2,const TestFunction &M, const TestFunction& N,double w,int DTM) const;
    
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  };
  
  class MyKoppelEquation : public virtual Equation
  {
  protected:

    mutable FemFunction* u1;
     mutable FemFunction* u2;
    
    double epsilon;
    mutable FemFunction* oldh;
  public:
    
    MyKoppelEquation() { abort(); }
    MyKoppelEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldh") != q.end() );
      oldh = &q["oldh"];
       assert(q.find("u1") != q.end() );
      u1= &q["u1"];
       assert(q.find("u2") != q.end() );
      u2= &q["u2"];
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    
    void Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const;
    void Nonlinear_Matrix(EntryMatrix&A, double s, const FemFunction &U1, const FemFunction& U2, const TestFunction &M, const TestFunction& N,double w,int DTM) const;
  
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  };
  
 class MyKoppelDualEquation : public virtual Equation
  {
  protected:

    mutable FemFunction* oldw;
    mutable FemFunction* dtu1;
    mutable FemFunction* dtu2;
    mutable FemFunction* dtu3;
    mutable FemFunction* u1;
    mutable FemFunction* u2;
    mutable FemFunction* u3;
    mutable FemFunction* oldz;
    mutable FemFunction* z;
    
    double epsilon;
  public:
    
    MyKoppelDualEquation() { abort(); }
    MyKoppelDualEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyDualTransortEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldw") != q.end() );
      oldw = &q["oldw"];
      assert(q.find("u2") != q.end() );
       u2 = &q["u2"];
        assert(q.find("dtu2") != q.end() );
       dtu2 = &q["dtu2"];
       
       
      if (!LASTDUAL){
        assert(q.find("u1") != q.end() );
       u1 = &q["u1"];  
      assert(q.find("dtu1") != q.end() );
      dtu1 = &q["dtu1"];
       assert(q.find("oldz") != q.end() );
      oldz= &q["oldz"];    
      
    }
      else
      {
      u1=NULL;    
	dtu1=NULL;
    oldz=NULL;
      } 
      
      if(!FIRSTDUAL){
          
    assert(q.find("u3") != q.end() );
      u3 = &q["u3"];
       assert(q.find("dtu3") != q.end() );
      dtu3 = &q["dtu3"];
       assert(q.find("z") != q.end() );
      z= &q["z"];   
      }
      else{
          
    u3=NULL;
    dtu3=NULL;
    z=NULL;   
    }
     
      
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Kopplung(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& OLDZ, const TestFunction& N,double w) const;
    void Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const;
    void Nonlinear_Matrix(EntryMatrix&A, double s, const FemFunction &U1, const FemFunction& U2,const TestFunction &M, const TestFunction& N,double w,int DTM) const;
    
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  };
  

  
 
}


/*----------------------------   equation.h     ---------------------------*/
/* end of #ifndef __equation_H */
#endif
/*----------------------------   equation.h     ---------------------------*/
