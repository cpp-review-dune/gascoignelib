#ifndef __numerusclaususadaptor
#define __numerusclaususadaptor

#include  "nvector.h"
#include  "paramfile.h"
#include  <string>

/*-----------------------------------------------------------*/

class NCAdaptor
{
protected:

  int    _n;
  double _p;
  double  etasum;
  const  nvector<double>&   eta;

public:

  NCAdaptor(const Gascoigne::ParamFile* paramfile, const nvector<double>& eta);
  void refine(nvector<int>& ref, nvector<int>& coarse) const;
};

#endif
