#ifndef __numerusclaususadaptor
#define __numerusclaususadaptor

#include  "nvector.h"
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

  NCAdaptor(const std::string& filename, const nvector<double>& eta);
  void refine(nvector<int>& ref, nvector<int>& coarse) const;
};

#endif
