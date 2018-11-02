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
  


namespace Gascoigne
{
  template class CGDiscQ12d;
  template class CGDiscQ22d;
  template class CGDiscQ13d;
  template class CGDiscQ23d;

  template class CGDiscQ23dLps;

}

