#include  "finiteelement.h"
#include  "finiteelement.xx"
#include  "transformation2d.h"
#include  "baseq12d.h"

/*-----------------------------------------------------*/

typedef Transformation2d<BaseQ12d>  TQ1_2D;

template class FiniteElement<2,1,TQ1_2D,BaseQ12d>;

/*-----------------------------------------------------*/
