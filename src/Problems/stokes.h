#ifndef  __Stokes_h
#define  __Stokes_h

#include  "equation.h"
#include  <string>

/*-----------------------------------------*/

class Stokes : public virtual Equation
{
protected:

  double visc;
  double penalty;

  double Laplace(const DerivativeVector& U, const TestFunction& N) const;
  double Divergence(const Gascoigne::FemFunction& U) const;

public:

  ~Stokes();
  Stokes();
  Stokes(const Gascoigne::ParamFile* pf);

  std::string GetName() const { return "Stokes";}

  int ncomp  () const { return 3; }
  //
  // Time
  //

  void SetTimePattern(TimePattern& P) const;
  
  void point(double _h, const Gascoigne::FemFunction& U, const Vertex2d& v) const {}

  //
  // Semilinear Form
  //

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};

#endif
