#ifndef  __StokesGls2d_h
#define  __StokesGls2d_h

#include  "stokes2d.h"
#include  "glsequation.h"
#include  "stabilization.h"

/*-----------------------------------------*/

class StokesGls2d : public Stokes2d, public virtual GlsEquation
{
protected:

  mutable Stabilization ST;

public:

  ~StokesGls2d();
  StokesGls2d();
  StokesGls2d(const Gascoigne::ParamFile* pf);

  std::string GetName() const { return "StokesGls2d";}

  //
  /// Computation of stabilization parameters
  //
  void glspoint(double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const;
  //
  /// for Galerkin-Least-Squares
  //
  void L(nvector<double>& dst, const Gascoigne::FemFunction& U) const;
  void S(nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;
};

#endif
