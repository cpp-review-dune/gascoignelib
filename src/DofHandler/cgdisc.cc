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

template<>
std::string
CGDiscQ12d::GetName() const
{
  return "Q12d";
}
template class CGDisc<2, 1, FiniteElementQ12d, ElementIntegratorQ12d>;

template<>
std::string
CGDiscQ22d::GetName() const
{
  return "Q22d";
}
template class CGDisc<2, 2, FiniteElementQ22d, ElementIntegratorQ22d>;

template<>
std::string
CGDiscQ13d::GetName() const
{
  return "Q13d";
}
template class CGDisc<3, 1, FiniteElementQ13d, ElementIntegratorQ13d>;

template<>
std::string
CGDiscQ23d::GetName() const
{
  return "Q23d";
}
template class CGDisc<3, 2, FiniteElementQ23d, ElementIntegratorQ23d>;

template class CGDisc<2, 1, FiniteElementQ22d, ElementIntegratorP12d>;

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

template class CGDisc<
  3,
  2,
  FiniteElement<3, 2, Transformation3d<BaseQ13dPatch>, BaseQ13dPatch>,
  ElementLpsIntegratorQ13d>;

template class CGDisc<3,
                      2,
                      FiniteElement<3, 2, Transformation3d<BaseQ23d>, BaseQ23d>,
                      ElementLpsIntegratorQ23d>;

template class CGDisc<
  2,
  2,
  FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>, BaseQ12dPatch>,
  ElementLpsIntegratorQ12d>;
template class CGDisc<2,
                      2,
                      FiniteElement<2, 1, Transformation2d<BaseQ22d>, BaseQ22d>,
                      ElementLpsIntegratorQ22d>;
template class CGDisc<
  2,
  2,
  FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>, BaseQ12dPatch>,
  ElementIntegratorQ12dPatch>;

} // namespace Gascoigne
