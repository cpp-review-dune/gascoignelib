//#include "cgdisc.xx"
#include "lagrangedisc.h"

#include "../Discretization/Q1/finiteelement.xx"

#include "../Common/compose_name.h"
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

#include "cgbase.h"
#include "elementintegrator.h"
#include "elementlpsintegrator.h"

namespace Gascoigne {
template class FiniteElement<2,
                             1,
                             Transformation2d<CGBase<2, 2>>,
                             CGBase<2, 2>>;
template class FiniteElement<2,
                             1,
                             Transformation2d<CGBase<2, 3>>,
                             CGBase<2, 3>>;
template class FiniteElement<2,
                             1,
                             Transformation2d<CGBase<2, 5>>,
                             CGBase<2, 5>>;
template class FiniteElement<3,
                             2,
                             Transformation3d<CGBase<3, 2>>,
                             CGBase<3, 2>>;
template class FiniteElement<3,
                             2,
                             Transformation3d<CGBase<3, 3>>,
                             CGBase<3, 3>>;
template class FiniteElement<3,
                             2,
                             Transformation3d<CGBase<3, 5>>,
                             CGBase<3, 5>>;

template class LagrangeDiscQ12d;
template class LagrangeDiscQ22d;
template class LagrangeDiscQ42d;

template class LagrangeDiscQ13d;
template class LagrangeDiscQ23d;
template class LagrangeDiscQ43d;
} // namespace Gascoigne
