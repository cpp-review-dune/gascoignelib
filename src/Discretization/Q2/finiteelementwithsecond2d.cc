#include  "transformation2d.h"
#include  "baseq22dwithsecond.h"
#include  "../Q1/finiteelement.xx"
#include  "finiteelementwithsecond.xx"
#include  "finiteelementwithsecond2d.xx"

/*-----------------------------------------------------*/

namespace Gascoigne
{
template class FiniteElement          <2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond>;
template class FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond>;
}
/*-----------------------------------------------------*/

