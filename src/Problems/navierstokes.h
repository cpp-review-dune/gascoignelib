#ifndef  __NavierStokes_h
#define  __NavierStokes_h

#include  "equation.h"

/*-----------------------------------------*/

class NavierStokes : public virtual Equation
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

  ~NavierStokes();
  NavierStokes();
  NavierStokes(const Gascoigne::ParamFile* pf);

  void OperatorStrong(Vector& b, const Gascoigne::FemFunction& U) const;

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
  void point(double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const { _h = h;}

  //
  // Semilinear Form
  //

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
};

#endif
