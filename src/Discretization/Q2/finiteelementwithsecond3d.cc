#include  "transformation3d.h"
#include  "baseq23dwithsecond.h"
#include  "finiteelement.xx"
#include  "finiteelementwithsecond.xx"
#include  "finiteelementwithsecond3d.xx"

/*-----------------------------------------------------*/

template class FiniteElement          <3, 2, Transformation3d<BaseQ23dWithSecond>, BaseQ23dWithSecond>;
template class FiniteElementWithSecond<3, 2, Transformation3d<BaseQ23dWithSecond>, BaseQ23dWithSecond>;

/*-----------------------------------------------------*/
