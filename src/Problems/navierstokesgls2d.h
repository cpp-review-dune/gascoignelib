#ifndef  __NavierStokesGls2d_h
#define  __NavierStokesGls2d_h

#include  "glsequation.h"
#include  "navierstokes2d.h"
#include  "glsstabilization.h"

/*-----------------------------------------*/

class NavierStokesGls2d : public NavierStokes2d, public virtual GlsEquation
{
protected:

  mutable GlsStabilization ST;

public:

  ~NavierStokesGls2d();
  NavierStokesGls2d();
  NavierStokesGls2d(const Gascoigne::ParamFile* pf);

  std::string GetName() const { return "NavierStokesGls2d";}

  void SetTime(double k)   { ST.DeltaT() = k;}
  //
  /// Computation of gls stabilization parameters
  //
  void glspoint(double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const;

  //
  /// for Galerkin-Least-Squares
  //
  void L(nvector<double>& dst, const Gascoigne::FemFunction& U) const;
  void S(nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void SMatrix(nvector<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::FemFunction& M, const Gascoigne::FemFunction& N) const;
};

#endif
