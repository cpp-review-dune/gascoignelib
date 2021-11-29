#include "cgdisc.xx"

#include "baseq12d.h"
#include "baseq13d.h"
#include "baseq22d.h"
#include "baseq23d.h"
#include "elementintegrator.h"
#include "elementlpsintegrator.h"
#include "finiteelement.h"
#include "integrationformula.h"
#include "transformation2d.h"
#include "transformation3d.h"

#include "baseq13dpatch.h"
#include "baseq1patch.h"

#include "finiteelement.xx"
#include "patchintegrationformula.h"

namespace Gascoigne {

template class CGDiscP12d;

template class CGDiscQ12d;
template class CGDiscQ22d;
template class CGDiscQ13d;
template class CGDiscQ23d;

template class FiniteElement<2,
                             1,
                             Transformation2d<BaseQ12dPatch>,
                             BaseQ12dPatch>;
template class FiniteElement<3,
                             2,
                             Transformation3d<BaseQ13dPatch>,
                             BaseQ13dPatch>;

template class FiniteElement<2, 1, Transformation2d<BaseQ22d>, BaseQ22d>;
template class FiniteElement<3, 2, Transformation3d<BaseQ23d>, BaseQ23d>;

template class CGDiscQ13dLps;
template class CGDiscQ23dLps;

template class CGDiscQ12dLps;
template class CGDiscQ22dLps;

} // namespace Gascoigne
