#ifndef  __Extrapolator_h
#define  __Extrapolator_h

#include  "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class Extrapolator
{
protected:

  std::vector<DoubleVector> vals;
  
  DoubleVector  valextra, order;

public:

  Extrapolator();
  ~Extrapolator();

  void Print();
  void NewValues(const DoubleVector& J);

};
}

/*-----------------------------------------*/

#endif
