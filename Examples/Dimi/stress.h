/*----------------------------   stress.h     ---------------------------*/
/*      $Id: myequation.h,v 1.1 2009/10/26 08:42:26 richter Exp $                 */
#ifndef __stress_H
#define __stress_H
/*----------------------------   stress.h     ---------------------------*/


#define  __MyEquation_h

#include  "paramfile.h"
#include  "equation.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class Stress : public virtual Equation
  {
  protected:
    mutable FemFunction* V;
    double rho;
     double  Tref, Lref,C, ellipse,Pstern;
    double deltamin;
    

  public:

    void SetFemData(FemData& q) const 
    {
      assert(q.find("v") != q.end() );
      V = &q["v"];
    }

    
    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    Stress() { abort(); }
    Stress(const ParamFile* pf);

    std::string GetName()  const { return "Stress";}

    int         GetNcomp() const {return 4;}

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
  };
}


/*----------------------------   equation.h     ---------------------------*/
/* end of #ifndef __equation_H */
#endif
/*----------------------------   equation.h     ---------------------------*/
