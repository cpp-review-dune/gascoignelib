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
#include  "paramfile.h"

class BackwardEquation : public Gascoigne::Equation
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


  std::string GetName() const { return "Local";}
  int  GetNcomp      () const {return 1;}
  void SetTimePattern(Gascoigne::TimePattern& P) const;
  
  void SetFemData(Gascoigne::FemData& Q) const;

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;

};


#endif
