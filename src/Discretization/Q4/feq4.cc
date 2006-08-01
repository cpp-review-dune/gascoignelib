#include "finiteelement.h"
#include "../Q1/finiteelement.xx"
#include "baseq12d.h"
#include "baseq13d.h"
#include "baseq22d.h"
#include "baseq23d.h"
#include "baseq42d.h"
#include "baseq43d.h"
#include "baseq2patch.h"
#include "baseq23dpatch.h"
#include "transformation2d.h"
#include "transformation3d.h"

namespace Gascoigne
{

/**********************************************************/

  typedef Transformation2d<BaseQ12d> TQ1_2D;
  typedef Transformation2d<BaseQ22d> TQ2_2D;
  typedef Transformation2d<BaseQ42d> TQ4_2D;

  template class FiniteElement<2,1,TQ1_2D,BaseQ42d>;
  template class FiniteElement<2,1,TQ2_2D,BaseQ42d>;
  template class FiniteElement<2,1,TQ4_2D,BaseQ42d>;
  template class FiniteElement<2,1,TQ1_2D,BaseQ22dPatch>;
  template class FiniteElement<2,1,TQ2_2D,BaseQ22dPatch>;

/**********************************************************/

  typedef Transformation3d<BaseQ13d> TQ1_3D;
  typedef Transformation3d<BaseQ23d> TQ2_3D;
  typedef Transformation3d<BaseQ43d> TQ4_3D;

  template class FiniteElement<3,2,TQ1_3D,BaseQ43d>;
  template class FiniteElement<3,2,TQ2_3D,BaseQ43d>;
  template class FiniteElement<3,2,TQ4_3D,BaseQ43d>;
  template class FiniteElement<3,2,TQ1_3D,BaseQ23dPatch>;
  template class FiniteElement<3,2,TQ2_3D,BaseQ23dPatch>;

/**********************************************************/

}

