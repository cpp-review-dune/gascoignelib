#ifndef  __NavierStokes2d_h
#define  __NavierStokes2d_h

#include  "equation.h"

/*-----------------------------------------*/

class NavierStokes2d : public virtual Equation
{
protected:

  mutable double _h, visc;
  double penalty,cut;

  double Laplace(const Gascoigne::TestFunction& U, 
		 const Gascoigne::TestFunction& N) const;
  
  double Convection(const Gascoigne::FemFunction& U, 
		    const Gascoigne::TestFunction& N) const;

  double Divergence(const Gascoigne::FemFunction& U) const;

public:

  ~NavierStokes2d();
  NavierStokes2d();
  NavierStokes2d(const Gascoigne::ParamFile* pf);

  void OperatorStrong(Gascoigne::DoubleVector& b, const Gascoigne::FemFunction& U) const;

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
  void point(double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const { _h = h;}

  //
  // Semilinear Form
  //

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
};

#endif
