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

public:


//
///  Constructor 
//
  LocalEquation(const Gascoigne::ParamFile* paramfile);

  std::string GetName() const { return "Local";}
  int    ncomp      () const {return 1;}
  void SetTimePattern(TimePattern& P) const;

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;

};


#endif
