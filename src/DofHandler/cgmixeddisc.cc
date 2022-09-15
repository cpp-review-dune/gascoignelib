#include "cgmixeddisc.h"

#include "baseq12d.h"
#include "baseq13d.h"
#include "baseq22d.h"
#include "baseq23d.h"
#include "mixedelementintegrator.h"

#include "finiteelement.h"
#include "integrationformula.h"
#include "transformation2d.h"
#include "transformation3d.h"

#include "baseq13dpatch.h"
#include "baseq1patch.h"

#include "finiteelement.xx"
#include "patchintegrationformula.h"

namespace Gascoigne {

template class Transformation2d<BaseQ12dPatch>;
//  template class FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>,
//  BaseQ22d>;

template class CGMixedDisc<
  2,
  2,
  FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>, BaseQ12dPatch>,
  FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>, BaseQ22d>,
  MixedElementIntegratorQ12dPatch>;
} // namespace Gascoigne
