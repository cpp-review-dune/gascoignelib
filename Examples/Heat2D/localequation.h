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
  LocalEquation(const Gascoigne::ParamFile* paramfile);

  std::string GetName() const { return "Local";}
  int    ncomp      () const {return 2;}
  void SetTimePattern(TimePattern& P) const;

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;
s  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;

  double GetUs() const {return _us;}
  double GetVs() const {return _vs;}

};


#endif
