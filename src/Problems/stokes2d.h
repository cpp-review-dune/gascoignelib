#ifndef  __Stokes2d_h
#define  __Stokes2d_h

#include  "equation.h"
#include  <string>

/*-----------------------------------------*/

namespace Gascoigne
{
class Stokes2d : public virtual Equation
{
protected:

  double visc;
  double penalty;

  double Laplace(const TestFunction& U, const TestFunction& N) const;
  double Divergence(const FemFunction& U) const;

public:

  ~Stokes2d();
  Stokes2d();
  Stokes2d(const ParamFile* pf);

  std::string GetName() const { return "Stokes2d";}

  int ncomp  () const { return 3; }
  //
  // Time
  //

  void SetTimePattern(TimePattern& P) const;
  
  void point(double _h, const FemFunction& U, const Vertex2d& v) const {}

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};
}

#endif
