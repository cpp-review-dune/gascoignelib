#ifndef __optadaptor_h
#define __optadaptor_h

#include "nvector.h"
#include "adaptordata.h"

/*********************************************************************/

class OptAdaptor
{
 protected:

  int d,p,p2,p4,n_aimed;
  int refined, double_refined, coarsened, marge, used;

  double co, dd, pp, factor;

  AdaptorData&             info;
  const nvector<double>&   vol;
        nvector<double>&   eta;

  void prepare();

public:

  OptAdaptor  (AdaptorData&, nvector<double>&, const nvector<double>&);

  void refine (nvector<int>&);
  void coarse (nvector<int>&);
  void RefineGnuplot (nvector<int>&);
};

/*********************************************************************/

#endif
