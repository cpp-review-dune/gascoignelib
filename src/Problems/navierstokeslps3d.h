#ifndef  __NavierStokesLps3d_h
#define  __NavierStokesLps3d_h

#include  "navierstokes3d.h"
#include  "lpsequation.h"
#include  "lpsstabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class NavierStokesLps3d : public NavierStokes3d, public virtual LpsEquation
{
  protected:

  mutable LpsStabilization ST;
  
  public:

  NavierStokesLps3d(const Gascoigne::ParamFile* filename);

  std::string GetName() const;

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of lps stabilization parameters
  //
  void lpspoint(double h, const FemFunction& U, const Vertex3d& v) const;
  //
  /// for local-projection stabilization (lps)
  //
  void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    
  void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;
};
}

#endif
