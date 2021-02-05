/*----------------------------   loop.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/

#include "boundaryfunction.h"
#include "gascoignemesh2d.h"
#include "local.h"
#include "meshagent.h"
#include "stdloop.h"
#include "usefullfunctionsbd.h"
#include "vertex.h"

#include "gascoignehash.h"
#include "multilevelsolver.h"
#include "solvers.h"
#include "stdmultilevelsolver.h"

namespace Gascoigne {

template <int DIM> class Loop : public StdLoop {

public:
  void BasicInit(const ParamFile *paramfile, const ProblemContainer *PC,
                 const FunctionalContainer *FC) {
    GetMeshAgentPointer() = new ProjectionOnFineMeshAgent(paramfile);
    // GetMultiLevelSolverPointer() = new StdMultiLevelSolver;
    GetMultiLevelSolverPointer() = new FSIMultiLevelSolver<DIM>;

    StdLoop::BasicInit(paramfile, PC, FC);
  }

  void run(const std::string &problemlabel);
};

} // namespace Gascoigne

/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
