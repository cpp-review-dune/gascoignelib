#ifndef  __StokesGls2d_h
#define  __StokesGls2d_h

#include  "stokes2d.h"
#include  "glsequation.h"
#include  "stabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class StokesGls2d : public Stokes2d, public virtual GlsEquation
{
protected:

  mutable Stabilization ST;

public:

  ~StokesGls2d();
  StokesGls2d();
  StokesGls2d(const ParamFile* pf);

  std::string GetName() const { return "StokesGls2d";}

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of stabilization parameters
  //
  void glspoint(double h, const FemFunction& U, const Vertex2d& v) const;
  //
  /// for Galerkin-Least-Squares
  //
  void L(DoubleVector& dst, const FemFunction& U) const;
  void S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;
};
}

#endif
