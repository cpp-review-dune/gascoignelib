#include  "finiteelement.h"
#include  "../Q1/finiteelement.xx"
#include  "transformation2d.h"
#include  "transformation3d.h"
#include  "baseq12d.h"
#include  "baseq22d.h"
#include  "baseq23d.h"
#include  "baseq1patch.h"
#include  "baseq13dpatch.h"

namespace Gascoigne
{
/*-----------------------------------------------------*/

typedef Transformation2d<BaseQ12d>  TQ1_2D;
typedef Transformation2d<BaseQ22d>  TQ2_2D;

template FiniteElement<2,1,TQ1_2D,BaseQ22d>;
template FiniteElement<2,1,TQ2_2D,BaseQ22d>;
template FiniteElement<2,1,TQ1_2D,BaseQ12dPatch>;

/*-----------------------------------------------------*/

typedef Transformation3d<BaseQ13d>  TQ1_3D;
typedef Transformation3d<BaseQ23d>  TQ2_3D;

template FiniteElement<3,2,TQ2_3D,BaseQ23d>;
template FiniteElement<3,2,TQ1_3D,BaseQ13dPatch>;
}
/*-----------------------------------------------------*/
