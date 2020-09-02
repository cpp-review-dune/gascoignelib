#include "cgdisc.xx"

#include "finiteelement.h"
#include "elementintegrator.h"
#include "elementlpsintegrator.h"
#include "integrationformula.h"
#include "transformation2d.h"
#include "transformation3d.h"
#include "baseq12d.h"
#include "baseq13d.h"
#include "baseq22d.h"
#include "baseq23d.h"

#include "baseq1patch.h"
#include "baseq13dpatch.h"

#include "patchintegrationformula.h"
#include "finiteelement.xx"

namespace Gascoigne
{
  template class CGDiscQ12d;
  template class CGDiscQ22d;
  template class CGDiscQ13d;
  template class CGDiscQ23d;

  template class FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>, BaseQ12dPatch>;
  template class FiniteElement<3, 2, Transformation3d<BaseQ13dPatch>, BaseQ13dPatch>;

  template class CGDiscQ13dLps;
  template class CGDiscQ23dLps;

  template class CGDiscQ12dLps;
  template class CGDiscQ22dLps;

}
