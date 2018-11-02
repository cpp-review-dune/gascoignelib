/*----------------------------   cglpsdisc.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __cglpsdisc_H
#define __cglpsdisc_H
/*----------------------------   cglpsdisc.h     ---------------------------*/

#include "cgdisc.h"

namespace Gascoigne
{
  
  // DIM=2,3
  // DEGREE = 1 (Q1) 2 (Q2)
  template <int DIM, int DEGREE, class FINITEELEMENT, class INTEGRATOR, class LPSINTEGRATOR>
  class CGLpsDisc : public CGDisc<DIM,DEGREE,FINITEELEMENT,INTEGRATOR,LPSINTEGRATOR>
  {
  };
}



/*----------------------------   cglpsdisc.h     ---------------------------*/
/* end of #ifndef __cglpsdisc_H */
#endif
/*----------------------------   cglpsdisc.h     ---------------------------*/
