#ifndef __optadaptor_h
#define __optadaptor_h

#include "gascoigne.h"
#include "adaptordata.h"

/*********************************************************************/

namespace Gascoigne
{
class OptAdaptor
{
 protected:

  int d,p,p2,p4,n_aimed;
  int refined, double_refined, coarsened, marge, used;

  double co, dd, pp, factor;

  AdaptorData&             info;
  const DoubleVector&   vol;
        DoubleVector&   eta;

  void prepare();

public:

  OptAdaptor  (AdaptorData&, DoubleVector&, const DoubleVector&);

  void refine (IntVector&);
  void coarse (IntVector&);
  void RefineGnuplot (IntVector&);
};
}

/*********************************************************************/

#endif
