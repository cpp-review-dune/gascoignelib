#ifndef  __StokesGls3d_h
#define  __StokesGls3d_h

#include  "stokes3d.h"
#include  "glsequation.h"
#include  "stabilization.h"

/*-----------------------------------------*/

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
  StokesGls3d(const Gascoigne::ParamFile* pf);

  std::string GetName() const { return "StokesGls3d";}

  //
  /// Computation of gls stabilization parameters
  //
  void glspoint(double h, const Gascoigne::FemFunction& U, const Vertex3d& v) const;
  //
  /// for Galerkin-Least-Squares
  //
  void L(nvector<double>& dst, const Gascoigne::FemFunction& U) const;
  void S(nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

};

#endif
