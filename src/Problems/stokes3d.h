#ifndef  __Stokes3d_h
#define  __Stokes3d_h

#include  "stokes.h"

/*-----------------------------------------*/

class Stokes3d : public Stokes
{
 protected:

  double Divergence(const Gascoigne::FemFunction& U) const;
  double Laplace(const DerivativeVector& U, const Gascoigne::TestFunction& N) const;

public:

  ~Stokes3d();
  Stokes3d();
  Stokes3d(const Gascoigne::ParamFile* pf);

  std::string GetName() const { return "Stokes3d";}

  int    ncomp  () const { return 4; }

  void point(double h, const Gascoigne::FemFunction& U, const Vertex3d& v) const {};

  //
  // Semilinear Form
  //

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
};

#endif
