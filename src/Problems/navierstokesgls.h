#ifndef  __NavierStokesGls_h
#define  __NavierStokesGls_h

#include  "glsequation.h"
#include  "navierstokes.h"
#include  "glsstabilization.h"

/*-----------------------------------------*/

class NavierStokesGls : public NavierStokes, public virtual GlsEquation
{
protected:

  mutable GlsStabilization ST;

public:

  ~NavierStokesGls();
  NavierStokesGls();
  NavierStokesGls(const std::string& paramfile);

  std::string GetName() const { return "NavierStokesGls";}

  void SetTime(double k)   { ST.DeltaT() = k;}
  //
  /// Computation of gls stabilization parameters
  //
  void glspoint(double h, const FemFunction& U, const Vertex2d& v) const;

  //
  /// for Galerkin-Least-Squares
  //
  void L(nvector<double>& dst, const FemFunction& U) const;
  void S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void SMatrix(nvector<double>& dst, const FemFunction& U, const FemFunction& M, const FemFunction& N) const;
};

#endif
