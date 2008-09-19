#ifndef  __NavierStokesSplit_h
#define  __NavierStokesSplit_h

#include  "lpsequation.h"
#include  "lpsstabilization.h"
#include  "paramfile.h"

/*---------------------------------------------------*/

class NavierStokesSplitLps2d : public virtual Gascoigne::LpsEquation
{
protected:

  mutable double _h, _visc;
  mutable Gascoigne::LpsStabilization ST;

  double Laplace(const Gascoigne::TestFunction& U, const Gascoigne::TestFunction& N) const;
  double Convection(const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;
  double Divergence(const Gascoigne::FemFunction& U) const;

public:

  NavierStokesSplitLps2d(const Gascoigne::ParamFile* filename);

  std::string GetName() const { return "NavierStokesSplitLps2d";}
  int  GetNcomp()       const { return 2; }

  void SetTime(double time, double dt) const;

  void SetTimePattern(Gascoigne::TimePattern& P) const;

  void point(double h, const Gascoigne::FemFunction& U, const Gascoigne::Vertex2d& v) const 
  { _h = h;}

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N) const;

  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
  //
  /// Computation of lps stabilization parameters
     //
    void lpspoint(double h, const Gascoigne::FemFunction& U, 
		  const Gascoigne::Vertex2d& v) const
  {
    ST.ReInit(h,_visc,U[0].m(),U[1].m());
  }
  //
  /// for local-projection stabilization (lps)
  //

  void StabForm(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
		const Gascoigne::FemFunction& UP, const Gascoigne::TestFunction& N) const;
    
  void StabMatrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
		  const Gascoigne::TestFunction& Np, const Gascoigne::TestFunction& Mp) const;
};

/*---------------------------------------------------*/

#endif
