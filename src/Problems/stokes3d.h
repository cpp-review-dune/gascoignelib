#ifndef  __Stokes3d_h
#define  __Stokes3d_h

#include  "stokes2d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class Stokes3d : public Stokes2d
{
 protected:

  double Divergence(const FemFunction& U) const;
  double Laplace(const TestFunction& U, const TestFunction& N) const;

public:

  ~Stokes3d();
  Stokes3d();
  Stokes3d(const ParamFile* pf);

  std::string GetName() const { return "Stokes3d";}

  int  GetNcomp  () const { return 4; }

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};
}

#endif
