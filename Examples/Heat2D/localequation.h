#ifndef  __LocalEquation_h
#define  __LocalEquation_h



/////////////////////////////////////////////
///
///@brief
///  ... comments LocalEquation

///
///
/////////////////////////////////////////////


#include  "equation.h"

class LocalEquation : public Equation
{
public:


private:


protected:

  mutable double _visc;
  mutable double _us, _vs;
  mutable double _h, _k, _r;

public:


//
///  Constructor 
//
  LocalEquation(const ParamFile* paramfile);

  std::string GetName() const { return "Local";}
  int    ncomp      () const {return 2;}
  void SetTimePattern(TimePattern& P) const;

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const DerivativeVector& M, const DerivativeVector& N) const;
  void Residual(Vector& b, const FemFunction& U, const DerivativeVector& N) const;


  double GetUs() const {return _us;}
  double GetVs() const {return _vs;}

};


#endif
