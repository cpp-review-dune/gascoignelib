#ifndef  __Functional_h
#define  __Functional_h

#include  <string>
#include  "equation.h"

//////////////////////////////////////////////
//
/// Interface class for Functional
//
//////////////////////////////////////////////

/*-----------------------------------------*/

class Functional : public Application
{
protected:

  typedef  nvector<double>         Vector;

  double  exact;
  bool    exactisknown;

  std::string beautifulname;

public:

  Functional() : Application(), exactisknown(0), exact(0.), beautifulname("NoBeautifulName")  {}
  ~Functional() {}
  Functional(const Functional& F) : Application(F)
    {
      exact = F.ExactValue();
      beautifulname = F.BeautifulName();
    } 

  virtual std::string GetName() const=0;

  double  ExactValue() const { return exact;}
  double& ExactValue()       { exactisknown = 1; return exact;}
  bool ExactValueIsKnown() const { return exactisknown; }

  std::string BeautifulName() const { return beautifulname;}
};


#endif
