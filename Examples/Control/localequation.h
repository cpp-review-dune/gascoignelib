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
  LocalEquation(const std::string& filename);

  std::string GetName() const { return "Local";}
  int    ncomp      () const {return 1;}
  void SetTimePattern(TimePattern& P) const;

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const DerivativeVector& M, const DerivativeVector& N) const;

};


#endif
