#include  "transformation2d.h"
#include  "baseq22dwithsecond.h"
#include  "finiteelement.xx"
#include  "finiteelementwithsecond.xx"
#include  "finiteelementwithsecond2d.xx"

/*-----------------------------------------------------*/

template class FiniteElement          <2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond>;
template class FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond>;

/*-----------------------------------------------------*/

