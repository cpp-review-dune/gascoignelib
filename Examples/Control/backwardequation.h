#ifndef  __BackwardEquation_h
#define  __BackwardEquation_h



/////////////////////////////////////////////
///
///@brief
///  ... comments BackwardEquation

///
///
/////////////////////////////////////////////


#include  "equation.h"

class BackwardEquation : public Equation
{
public:


private:


  mutable const Gascoigne::FemFunction* q;

protected:

  mutable double _visc;

public:


//
///  Constructor 
//
  BackwardEquation(const Gascoigne::ParamFile* paramfile);

  void point(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& Q, const Vertex2d& v) const;

  std::string GetName() const { return "Local";}
  int    ncomp      () const {return 1;}
  void SetTimePattern(TimePattern& P) const;

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const DerivativeVector& M, const DerivativeVector& N) const;

};


#endif
