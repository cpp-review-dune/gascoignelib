#ifndef  __Stokes2d_h
#define  __Stokes2d_h

#include  "equation.h"
#include  <string>

/*-----------------------------------------*/

class Stokes2d : public virtual Equation
{
protected:

  double visc;
  double penalty;

  double Laplace(const Gascoigne::TestFunction& U, const Gascoigne::TestFunction& N) const;
  double Divergence(const Gascoigne::FemFunction& U) const;

public:

  ~Stokes2d();
  Stokes2d();
  Stokes2d(const Gascoigne::ParamFile* pf);

  std::string GetName() const { return "Stokes2d";}

  int ncomp  () const { return 3; }
  //
  // Time
  //

  void SetTimePattern(TimePattern& P) const;
  
  void point(double _h, const Gascoigne::FemFunction& U, const Vertex2d& v) const {}

  //
  // Semilinear Form
  //

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
};

#endif