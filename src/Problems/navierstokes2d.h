#ifndef  __NavierStokes2d_h
#define  __NavierStokes2d_h

#include  "equation.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class NavierStokes2d : public virtual Equation
{
protected:

  mutable double _h, visc;
  double penalty,cut;

  double Laplace(const TestFunction& U, 
		 const TestFunction& N) const;
  
  double Convection(const FemFunction& U, 
		    const TestFunction& N) const;

  double Divergence(const FemFunction& U) const;

public:

  ~NavierStokes2d();
  NavierStokes2d();
  NavierStokes2d(const ParamFile* pf);

  void OperatorStrong(DoubleVector& b, const FemFunction& U) const;

  std::string GetName() const { return "NavierStokes2d";}
  double Getpenalty()const {return penalty;}

  int  ncomp() const { return 3; }

  //
  // Time
  //

  void SetTimePattern(TimePattern& P) const;

  //
  /// Computation of coefficients at each integration point.
  /// In the case of Navier-Stokes, this function is empty.
  //
  void point(double h, const FemFunction& U, const Vertex2d& v) const { _h = h;}

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};
}

#endif
