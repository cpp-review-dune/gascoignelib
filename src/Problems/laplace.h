#ifndef  __Laplace_h
#define  __Laplace_h

#include  "equation.h"

/*-----------------------------------------*/

class Laplace : public virtual Equation
{
  protected:
  
  double gamma;
  mutable double visc;
  
  public:

  Laplace();
  Laplace(const std::string& filename);


  std::string GetName() const { return "Laplace";}

  int    ncomp      () const {return 1;}

  void OperatorStrong(Vector& b, const FemFunction& U) const;

  void SetTimePattern(TimePattern& P) const;

  void point(double h, const FemFunction& U, const Vertex2d& v) const {}

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const DerivativeVector& M, const TestFunction& N) const;
};


#endif
