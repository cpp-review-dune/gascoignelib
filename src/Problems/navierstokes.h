#ifndef  __NavierStokes_h
#define  __NavierStokes_h

#include  "equation.h"

using namespace Gascoigne;

/*-----------------------------------------*/

class NavierStokes : public virtual Equation
{
protected:

  mutable double _h, visc;
  double penalty,cut;

  double Laplace(const DerivativeVector& U, 
		 const TestFunction& N) const;
  
  double Convection(const FemFunction& U, 
		    const TestFunction& N) const;

  double Divergence(const FemFunction& U) const;

public:

  ~NavierStokes();
  NavierStokes();
  NavierStokes(const ParamFile* pf);

  void OperatorStrong(Vector& b, const FemFunction& U) const;

  std::string GetName() const { return "NavierStokes";}
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

#endif
