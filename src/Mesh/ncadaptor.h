#ifndef __numerusclaususadaptor
#define __numerusclaususadaptor

#include  "nvector.h"
#include  "paramfile.h"
#include  <string>

using namespace Gascoigne;

/*-----------------------------------------------------------*/

class NCAdaptor
{
protected:

  int    _n;
  double _p;
  double  etasum;
  const  nvector<double>&   eta;

public:

  NCAdaptor(const ParamFile* paramfile, const nvector<double>& eta);
  void refine(nvector<int>& ref, nvector<int>& coarse) const;
};

#endif
