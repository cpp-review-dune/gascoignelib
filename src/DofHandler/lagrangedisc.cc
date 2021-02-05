//#include "cgdisc.xx"
#include "lagrangedisc.h"

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

#include "cgbase.h"
#include "compose_name.h"

namespace Gascoigne {
template class FiniteElement<2, 1, Transformation2d<CGBase<2, 2>>,
                             CGBase<2, 2>>;
template class FiniteElement<2, 1, Transformation2d<CGBase<2, 3>>,
                             CGBase<2, 3>>;
template class FiniteElement<2, 1, Transformation2d<CGBase<2, 5>>,
                             CGBase<2, 5>>;
template class FiniteElement<3, 2, Transformation3d<CGBase<3, 2>>,
                             CGBase<3, 2>>;
template class FiniteElement<3, 2, Transformation3d<CGBase<3, 3>>,
                             CGBase<3, 3>>;
template class FiniteElement<3, 2, Transformation3d<CGBase<3, 5>>,
                             CGBase<3, 5>>;

template class LagrangeDiscQ12d;
template class LagrangeDiscQ22d;
template class LagrangeDiscQ42d;

template class LagrangeDiscQ13d;
template class LagrangeDiscQ23d;
template class LagrangeDiscQ43d;
} // namespace Gascoigne
