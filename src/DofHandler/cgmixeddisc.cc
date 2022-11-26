#include "../Discretization/Q1/finiteelement.xx"

#include "cgmixeddisc.h"

#include "../Discretization/Q1/baseq12d.h"
#include "../Discretization/Q1/baseq13d.h"
#include "../Discretization/Q1/finiteelement.h"
#include "../Discretization/Q1/integrationformula.h"
#include "../Discretization/Q1/transformation2d.h"
#include "../Discretization/Q1/transformation3d.h"
#include "../Discretization/Q2/baseq13dpatch.h"
#include "../Discretization/Q2/baseq1patch.h"
#include "../Discretization/Q2/baseq22d.h"
#include "../Discretization/Q2/baseq23d.h"
#include "../Discretization/Q2/patchintegrationformula.h"

#include "mixedelementintegrator.h"

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
