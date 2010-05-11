#ifndef  __waveequation_h
#define  __waveequation_h

#include  "equation.h"
#include  "boundaryequation.h"

/*---------------------------------------------------*/

class WaveEquation : public Gascoigne::Equation
{
  double c2;
  double Laplace(const Gascoigne::TestFunction& U, const Gascoigne::TestFunction& N) const;

public:

  WaveEquation();

  int  GetNcomp() const { return 1; }

  std::string GetName() const { return "WaveEquation";}

  void SetTimePattern(Gascoigne::TimePattern& TP) const 
  { TP(0,0) = 1.;}

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N) const;

  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M, 
	      const Gascoigne::TestFunction& N) const;
};

/*---------------------------------------------------*/

class WaveBoundaryEquation : public Gascoigne::BoundaryEquation
{
  // for absorbing boundary condition
  //
  double c2;

public:

  WaveBoundaryEquation();
  int GetNcomp()const { return 1;}

  std::string GetName() const { return "WaveBoundaryEquation";}

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N, int color) const;

  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M,
              const Gascoigne::TestFunction& N, int color) const;
};

/*---------------------------------------------------*/

#endif
