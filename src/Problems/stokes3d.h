#ifndef  __Stokes3d_h
#define  __Stokes3d_h

#include  "stokes.h"

/*-----------------------------------------*/

class Stokes3d : public Stokes
{
 protected:

  double Divergence(const FemFunction& U) const;
  double Laplace(const DerivativeVector& U, const TestFunction& N) const;

public:

  ~Stokes3d();
  Stokes3d();
  Stokes3d(const std::string& filename);

  std::string GetName() const { return "Stokes3d";}

  int    ncomp  () const { return 4; }

  void point(double h, const FemFunction& U, const Vertex3d& v) const {};

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};

#endif
