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

extern bool Jump;
namespace Gascoigne
{
  
  class TransportEquation : public virtual Equation
  {
  protected:

    mutable FemFunction* oldH;
    mutable FemFunction* V;


   //     mutable FemFunction V;
    mutable double h_;
  public:
    
    TransportEquation() { abort(); }
    TransportEquation(const ParamFile* pf);

    std::string GetName()  const { return "TransportEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldH") != q.end() );
      oldH = &q["oldH"];
 
      assert(q.find("V") != q.end() );
      V = &q["V"];
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  };


  class MyDualEquation : public virtual Equation
  {
  protected:
    mutable FemFunction* nextV;
    mutable FemFunction* V;
    mutable FemFunction* nextQ;
    mutable FemFunction* nextZ;
    
  mutable FemFunction* oldH;

    double _split, f;
    
      mutable double uwx;
    mutable double uwy;
        double alpha0,  tau0, Cstern, Pstern, Tref, Lref, MZ, ellipse, C, Cdw, rhow, theta_w;
         mutable double rho,h_;
  public:
    
    MyDualEquation() { abort(); }
    MyDualEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyDualEquation";}
    int         GetNcomp() const {return 2;}
    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("nextZ") != q.end() ); 
      nextZ = &q["nextZ"];
      
      assert(q.find("nextV") != q.end() ); 
      nextV = &q["nextV"];

      assert(q.find("V") != q.end() );
      V= &q["V"];
      

      assert(q.find("oldH") != q.end() );
      oldH= &q["oldH"];
      
      
      assert(q.find("nextQ") != q.end() );
      nextQ= &q["nextQ"];
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& Z, const TestFunction& M, const TestFunction& N) const;

  };
  
   class MyBurgerEquation : public virtual Equation
  {
  protected:
      
    mutable FemFunction* oldV;
    mutable FemFunction* oldH;
  //   mutable FemFunction H;
    double epsilon, _split;
    mutable double h_;
    double alpha0,  tau0, Cstern, Pstern, Tref, Lref, MZ, ellipse, C, Cdw, rhow, theta_w,f;
    double deltamin, shock0;
    mutable double shock;
    mutable double alpha,tau;
    mutable double visc;
    mutable double uwx;
    mutable double uwy;
    
     mutable double rho;
    mutable double gamma;
    double vin0;
   
    mutable double v_in, gamma0;
    mutable Vertex2d _n;
  public:
    
    MyBurgerEquation() { abort(); }
    MyBurgerEquation(const ParamFile* pf);

    std::string GetName()  const { return "TransportEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldV") != q.end() );
      oldV = &q["oldV"];
      assert(q.find("oldH") != q.end() );
      oldH= &q["oldH"];
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
   
  };
  
 class MyBurgerDualEquation : public virtual Equation
  {
  protected:

     mutable FemFunction* nextQ;
     mutable FemFunction* V;
     mutable FemFunction* H;
     mutable FemFunction* oldH;
     mutable FemFunction* Z;
    //   mutable FemFunction H;
     
     double epsilon, _split;
    mutable double rho;
    double alpha0,  tau0, Cstern, Pstern, Tref, Lref, MZ, ellipse, C, Cdw, rhow, theta_w,f;
    double deltamin;
    mutable double uwx;
     mutable double uwy;
     
           mutable double h_;
       
     
  public:
    
    MyBurgerDualEquation() { abort(); }
    MyBurgerDualEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyDualEquation";}
    int         GetNcomp() const {return 2;}


    
    void SetFemData(FemData& q) const 
    {
      assert(q.find("V") != q.end() ); 
      V = &q["V"];
      
      assert(q.find("nextQ") != q.end() ); 
      nextQ = &q["nextQ"];
      
      assert(q.find("H") != q.end() ); 
      H = &q["H"];
      
      assert(q.find("oldH") != q.end() ); 
      oldH = &q["oldH"];

      assert(q.find("Z") != q.end() ); 
      Z = &q["Z"]; 
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& Z, const TestFunction& M, const TestFunction& N) const;

  };
  
  

}


/*----------------------------   equation.h     ---------------------------*/
/* end of #ifndef __equation_H */
#endif
/*----------------------------   equation.h     ---------------------------*/
