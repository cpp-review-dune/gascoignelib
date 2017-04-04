/*----------------------------   myequation.h     ---------------------------*/
/*      $Id: myequation.h,v 1.1 2012/09/02 16:43:19 richter Exp $                 */
#ifndef __myequation_H
#define __myequation_H
/*----------------------------   myequation.h     ---------------------------*/

#include  "equation.h"
#include  "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class MyEquation : public virtual Equation
  {
  protected:
  
    double __eps;
    mutable double __b0,__b1;
    
    
  public:

    MyEquation() { abort(); }
  
    MyEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyEquation";}

    int         GetNcomp() const {return 1;}

    void point(double h, const FemFunction& U, const Vertex2d& v) const;
  
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

    void SetTimePattern(TimePattern& P) const
    {
      P(0,0) = 1.0;
    }
    

  };
}



/*----------------------------   myequation.h     ---------------------------*/
/* end of #ifndef __myequation_H */
#endif
/*----------------------------   myequation.h     ---------------------------*/
