#ifndef  __NavierStokesGls2d_h
#define  __NavierStokesGls2d_h

#include  "glsequation.h"
#include  "navierstokes2d.h"
#include  "glsstabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class NavierStokesGls2d : public NavierStokes2d, public virtual GlsEquation
{
protected:

  mutable GlsStabilization ST;

public:

  ~NavierStokesGls2d();
  NavierStokesGls2d();
  NavierStokesGls2d(const ParamFile* pf);

  std::string GetName() const { return "NavierStokesGls2d";}

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of gls stabilization parameters
  //
  void glspoint(double h, const FemFunction& U, const Vertex2d& v) const;

  //
  /// for Galerkin-Least-Squares
  //
  void L(DoubleVector& dst, const FemFunction& U) const;
  void S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void SMatrix(DoubleVector& dst, const FemFunction& U, const FemFunction& M, const FemFunction& N) const;
};
}

#endif
