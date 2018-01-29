/*----------------------------   equation.h     ---------------------------*/
/*      $Id: myequation.h,v 1.1 2009/10/26 08:42:26 richter Exp $                 */
#ifndef __water_H
#define __water_H
/*----------------------------   equation.h     ---------------------------*/


#define  __MyEquation_h

#include  "paramfile.h"
#include  "lpsequation.h"
#include  "boundaryequation.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class Water : public virtual Equation, public BoundaryEquation
  {
  protected:
    double alpha0, visc,delta0;

    double Tref,Lref,vin0;
    
    mutable double v_in;
    mutable Vertex2d _n;

    // Parameter fuer die schwachen Randwerte gamma = gamma0/h
    double gamma0, lps0;
    mutable double gamma,lps;

    mutable double alpha;
    mutable double delta;

  public:

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    Water() { abort(); }
    Water(const ParamFile* pf);

    std::string GetName()  const { return "Water";}

    int         GetNcomp() const {return 3;}

    // Gleichung

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

    // Rand
    void pointboundary(double h, const FemFunction& U, const Vertex2d& v, const Vertex2d& n) const;
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const;
    void Matrix(EntryMatrix& E, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const;

    /* // Stabilisierung */
    /* void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const; */
    /* void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const; */
    /* void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const; */

    
  };
}


/*----------------------------   equation.h     ---------------------------*/
/* end of #ifndef __equation_H */
#endif
/*----------------------------   equation.h     ---------------------------*/
