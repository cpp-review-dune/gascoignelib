#ifndef  __StokesGls_h
#define  __StokesGls_h

#include  "stokes.h"
#include  "glsequation.h"
#include  "stabilization.h"

/*-----------------------------------------*/

class StokesGls : public Stokes, public virtual GlsEquation
{
protected:

  mutable Stabilization ST;

public:

  ~StokesGls();
  StokesGls();
  StokesGls(const ParamFile* pf);

  std::string GetName() const { return "StokesGls";}

  //
  /// Computation of stabilization parameters
  //
  void glspoint(double h, const FemFunction& U, const Vertex2d& v) const;
  //
  /// for Galerkin-Least-Squares
  //
  void L(nvector<double>& dst, const FemFunction& U) const;
  void S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;
};

#endif
