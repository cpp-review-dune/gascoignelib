#ifndef __numerusclaususadaptor
#define __numerusclaususadaptor

#include  "gascoigne.h"
#include  "paramfile.h"
#include  <string>

/*-----------------------------------------------------------*/

namespace Gascoigne
{
class NCAdaptor
{
protected:

  int    _n;
  double _p;
  double  etasum;
  const  DoubleVector&   eta;

public:

  NCAdaptor(const ParamFile* paramfile, const DoubleVector& eta);
  void refine(IntVector& ref, IntVector& coarse) const;
};
}

#endif
