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

namespace Gascoigne
{
  
  class MyEquation : public virtual Equation
  {
  protected:

    mutable FemFunction* oldu;
    mutable double beta_x,beta_y;
    double epsilon;
  public:
    
    MyEquation() { abort(); }
    MyEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyEquation";}
    int         GetNcomp() const {return 1;}


    
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
     mutable FemFunction* Pu;
     mutable double beta_x,beta_y;
     double epsilon;
  public:
    
    MyDualEquation() { abort(); }
    MyDualEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyDualEquation";}
    int         GetNcomp() const {return 1;}


    
    void SetFemData(FemData& q) const 
    {
       assert(q.find("oldz") != q.end() ); 
       oldz = &q["oldz"];
       
     assert(q.find("Pu") != q.end() ); 
     Pu = &q["Pu"];
       
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  };

class MyDualHelpEquation : public virtual Equation
  {
  protected:

    mutable FemFunction *U_h;
    mutable FemFunction *Uold;
    mutable double beta_x,beta_y;
    double epsilon;
    
  public:
    
    MyDualHelpEquation() { abort(); }
    MyDualHelpEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyDualHelpEquation";}
    int         GetNcomp() const {return 1;}

    void SetFemData(FemData& q) const 
    {
      assert(q.find("U_h") != q.end() );
      U_h = &q["U_h"];

      
      assert(q.find("Uold") != q.end() );
      Uold= &q["Uold"];
     
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
