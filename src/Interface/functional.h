#ifndef  __Functional_h
#define  __Functional_h

#include  <string>
#include  "derivativevector.h"
#include  "equation.h"

//////////////////////////////////////////////
//
/// Interface class for Functional
//
//////////////////////////////////////////////


/*-----------------------------------------*/


class Functional
{
protected:

  typedef  nvector<double>         Vector;
  typedef  DerivativeVector  TestFunction;

  double  exact;
  bool    exactisknown;

public:

  Functional() : exactisknown(0), exact(0.)  {}
  virtual ~Functional() {}

  virtual std::string GetName() const=0;

  double  ExactValue() const { return exact;}
  double& ExactValue()       { exactisknown = 1; return exact;}
  bool ExactValueIsKnown() const { return exactisknown; }
};


#endif
