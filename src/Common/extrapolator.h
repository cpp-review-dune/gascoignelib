#ifndef  __Extrapolator_h
#define  __Extrapolator_h

#include  "nvector.h"

/*-----------------------------------------*/

class Extrapolator
{
protected:

  typedef nvector<double> dvector;

  std::vector<dvector> vals;
  
  dvector  valextra, order;

public:

  Extrapolator();
  ~Extrapolator();

  void Print();
  void NewValues(const dvector& J);

};

/*-----------------------------------------*/

#endif
