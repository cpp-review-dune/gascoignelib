/*----------------------------   loop.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/

#include "stdloop.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include "vertex.h"
#include  "gascoignemesh2d.h"
#include  "solvers.h"



namespace Gascoigne
{

  template<int DIM>
    class Loop : public StdLoop
    {
      
    public:
      
      void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,const FunctionalContainer* FC)
      {
	GetMultiLevelSolverPointer() = new FSIMultiLevelSolver<DIM>;
	StdLoop::BasicInit(paramfile,PC,FC);
      }
    };
  
}


/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
