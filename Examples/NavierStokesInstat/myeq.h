/*----------------------------   myeq.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __myeq_H
#define __myeq_H
/*----------------------------   myeq.h     ---------------------------*/



#include  "equation.h"
#include  "lpsequation.h"
#include  "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class MyEQ : public virtual Equation, public LpsEquation
{
protected:

  mutable double _visc, _lps;
  double _lps0;

  mutable FemFunction* old;
  void SetFemData(FemData& q) const 
  { assert(q.find("old")!=q.end()); old = &q["old"]; }
  


  double Laplace(const TestFunction& U, const TestFunction& N) const;
  double Convection(const FemFunction& U, const TestFunction& N) const;

  double Divergence(const FemFunction& U) const;


  

public:

  ~MyEQ();
  MyEQ();
  MyEQ(const ParamFile* pf);

  void OperatorStrong(DoubleVector& b, const FemFunction& U) const;

  std::string GetName() const { return "MyEQ";}

  int  GetNcomp() const { return 3; }

  //
  // Time
  //

  void SetTimePattern(TimePattern& P) const;

  //
  /// Computation of coefficients at each integration point.
  /// In the case of Navier-Stokes, this function is empty.
  //
  void point(double h, const FemFunction& U, const Vertex2d& v) const { }

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
  
  //
  void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const
  {
    _lps = _lps0 * h * h/_visc;
  }
  void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
  {
    b[0] += _lps * (UP[0].x()*N.x()+UP[0].y()*N.y());
  }
  
    
  void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
  {
    A(0,0) += _lps * (Mp.x()*Np.x()+Mp.y()*Np.y());
  }
  
};

  
}



/*----------------------------   myeq.h     ---------------------------*/
/* end of #ifndef __myeq_H */
#endif
/*----------------------------   myeq.h     ---------------------------*/
