#include  "finiteelement.h"
#include  "transformation3d.h"
#include  "baseq13d.h"
#include  "finiteelement.xx"

/*-----------------------------------------------------*/

namespace Gascoigne
{
template class FiniteElement<3,2,Transformation3d<BaseQ13d>,BaseQ13d>;
}
