#ifndef  __NavierStokesLps2d_h
#define  __NavierStokesLps2d_h

#include  "lpsequation.h"
#include  "navierstokes2d.h"
#include  "lpsstabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class NavierStokesLps2d : public NavierStokes2d, public virtual LpsEquation
{
protected:

  mutable LpsStabilization ST;

public:

  NavierStokesLps2d(const Gascoigne::ParamFile* filename);

  std::string GetName() const { return "NavierStokesLps2d";}

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of lps stabilization parameters
  //
  void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const;
  //
  /// for local-projection stabilization (lps)
  //

  void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    
  void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;
};
}

#endif
