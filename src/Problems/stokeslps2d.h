#ifndef  __StokesLps2d_h
#define  __StokesLps2d_h

#include  "stokes2d.h"
#include  "lpsequation.h"
#include  "stabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class StokesLps2d : public Stokes2d, public virtual LpsEquation
{
protected:

  mutable Stabilization ST;

public:

  ~StokesLps2d();
  StokesLps2d();
  StokesLps2d(const ParamFile* filename);

  std::string GetName() const { return "StokesLps2d";}
  void SetTime(double k)   { ST.DeltaT() = k;}

  //
  /// Computation of stabilization parameters
  //
  void lpspoint(double h, const FemFunction& U, FemData& Q, const Vertex2d& v) const;

  void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    
  void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const DerivativeVector& Mp) const;

};
}

/*-----------------------------------------*/

#endif
