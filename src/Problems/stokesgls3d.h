#ifndef  __StokesGls3d_h
#define  __StokesGls3d_h

#include  "stokes3d.h"
#include  "glsequation.h"
#include  "stabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class StokesGls3d : public Stokes3d, public virtual GlsEquation
{
protected:

  //
  /// handles the stabilization parameters
  //
  mutable Stabilization ST;

public:

  ~StokesGls3d();
  StokesGls3d();
  StokesGls3d(const ParamFile* pf);

  std::string GetName() const { return "StokesGls3d";}

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of gls stabilization parameters
  //
  void glspoint(double h, const FemFunction& U, const Vertex3d& v) const;
  //
  /// for Galerkin-Least-Squares
  //
  void L(DoubleVector& dst, const FemFunction& U) const;
  void S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

};
}

#endif
