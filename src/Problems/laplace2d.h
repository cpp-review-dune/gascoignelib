#ifndef  __Laplace2d_h
#define  __Laplace2d_h

#include  "equation.h"

/*-----------------------------------------*/

class Laplace2d : public virtual Equation
{
  protected:
  
  double gamma;
  mutable double visc;
  
  public:

  Laplace2d();
  Laplace2d(const Gascoigne::ParamFile* pf);

  std::string GetName() const { return "Laplace2d";}

  int    ncomp      () const {return 1;}

  void OperatorStrong(Vector& b, const Gascoigne::FemFunction& U) const;

  void SetTimePattern(TimePattern& P) const;

  void point(double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const {}

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
};


#endif
