/*----------------------------   equation.h     ---------------------------*/
/*      $Id: myequation.h,v 1.1 2009/10/26 08:42:26 richter Exp $                 */
#ifndef __equation_H
#define __equation_H
/*----------------------------   equation.h     ---------------------------*/


#define  __MyEquation_h

#include  "paramfile.h"
#include  "equation.h"
#include "filescanner.h"

using namespace std;

/*-----------------------------------------*/

extern bool Jump;
namespace Gascoigne
{

  class Ice
  {
  public:
    
    double Pstern, Tref, Lref, MZ, ellipse, C, Cdw, rhow, theta_w,f;
    double deltamin;
    mutable double uwx;
    mutable double uwy;
    mutable double rho;

    Ice(){assert(0);}
    Ice(const ParamFile* pf)
      {
	DataFormatHandler DFH;
	DFH.insert("rho", &rho, 0.);
	DFH.insert("rhow", &rhow, 0.);
	DFH.insert("Tref",&Tref,0.0);
	DFH.insert("Lref",&Lref,0.0);
	DFH.insert("Pstern",&Pstern,2.75e4);
	DFH.insert("ellipse",&ellipse,2.0);
	DFH.insert("C",&C,20.0);
	DFH.insert("Cdw",&Cdw,5.2e-3);
	DFH.insert("f",&f,0.0);
	DFH.insert("theta_w",&theta_w,0.0);
	
	FileScanner FS(DFH);
	FS.NoComplain();
	FS.readfile(pf, "Equation");
	assert(rho>0);
	assert(Tref>0);
	assert(Lref>0);
	assert(Pstern>0);
	assert(ellipse>0);
	assert(C>0);
	assert(f>0);
	
	MZ = 0.5*Tref*Tref * Pstern / rho / Lref / Lref;
      }

    double delta(const FemFunction& U) const;
    double delta_0(const FemFunction& V, const TestFunction& M) const;
    double delta_1(const FemFunction& V, const TestFunction& M) const;
    double SIGMA_ij(const FemFunction& V, const TestFunction& M, const TestFunction& N, int i, int j) const;
  };
  

  
  class TransportEquation : public virtual Equation
  {
  protected:

    mutable FemFunction* oldH;
    mutable FemFunction* V;


    //     mutable FemFunction V;
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


  class MyDualEquation : public virtual Equation, public Ice
  {
  protected:
    mutable FemFunction* nextV;
    mutable FemFunction* V;
    mutable FemFunction* nextQ;
    mutable FemFunction* nextZ;
    
    mutable FemFunction* H;

    
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
      

      assert(q.find("H") != q.end() );
      H= &q["H"];
      
      
      assert(q.find("nextQ") != q.end() );
      nextQ= &q["nextQ"];
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& Z, const TestFunction& M, const TestFunction& N) const;

  };

  class MyBurgerEquation : public virtual Equation, public Ice
  {
  protected:
      
    mutable FemFunction* oldV;
    mutable FemFunction* oldH;
    //   mutable FemFunction H;


  
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
  
  class MyBurgerDualEquation : public virtual Equation, public Ice
  {
  protected:

    mutable FemFunction* nextQ;
    mutable FemFunction* V;
    mutable FemFunction* H;
    mutable FemFunction* oldH;
    mutable FemFunction* Z;
    //   mutable FemFunction H;
   
     
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
