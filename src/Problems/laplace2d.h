#ifndef  __Laplace2d_h
#define  __Laplace2d_h

#include  "equation.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class Laplace2d : public virtual Equation
{
  protected:
  
  double gamma;
  mutable double visc;
  
  public:

  Laplace2d();
  Laplace2d(const ParamFile* pf);

  std::string GetName() const { return "Laplace2d";}

  int    ncomp      () const {return 1;}

  void OperatorStrong(DoubleVector& b, const FemFunction& U) const;

  void SetTimePattern(TimePattern& P) const;

  void point(double h, const FemFunction& U, const Vertex2d& v) const {}

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};
}

#endif
