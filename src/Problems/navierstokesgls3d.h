#ifndef  __NavierStokesGls3d_h
#define  __NavierStokesGls3d_h

#include  "navierstokes3d.h"
#include  "glsequation.h"
#include  "glsstabilization.h"

/*-----------------------------------------*/

class NavierStokesGls3d : public NavierStokes3d, public virtual GlsEquation
{
  protected:

  mutable GlsStabilization ST;
  
  public:

  ~NavierStokesGls3d();
  NavierStokesGls3d();
  NavierStokesGls3d(const std::string& parmfile);

  std::string GetName() const;
  //
  /// Computation of gls stabilization parameters
  //
  void glspoint(double h, const FemFunction& U, const Vertex3d& v) const;
  void glspointmatrix(double _h, const FemFunction& U, FemData& Q, const Vertex3d& v) const;

  //
  /// for Galerkin-Least-Squares
  //
  void L(nvector<double>& dst, const FemFunction& U) const;
  void S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;
};

#endif
