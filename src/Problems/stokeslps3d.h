#ifndef  __StokesLps3d_h
#define  __StokesLps3d_h

#include  "stokes3d.h"
#include  "lpsequation.h"
#include  "stabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class StokesLps3d : public Stokes3d, public virtual LpsEquation
{
protected:

  //
  /// handles the stabilization parameters
  //
  mutable Stabilization ST;

public:

  ~StokesLps3d();
  StokesLps3d();
  StokesLps3d(const ParamFile* filename);

  std::string GetName() const { return "StokesLps3d";}

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of lps stabilization parameters
  //
  void lpspoint(double h, const FemFunction& U, const Vertex3d& v) const;

  //
  /// Stablization terms
  //
  void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    
  void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const DerivativeVector& Mp) const;
};
}

#endif
