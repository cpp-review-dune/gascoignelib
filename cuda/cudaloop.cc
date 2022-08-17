#include "cudaloop.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

#include "cudamultilevelsolver.h"

#include "dataformathandler.h"
#include "filescanner.h"
#include "nlinfo.h"
#include "solverinfos.h"
#include "stdmultilevelsolver.h"
#include "stopwatch.h"
#include "vectorinterface.h"

namespace Gascoigne {
class FunctionalContainer;
class ParamFile;
class ProblemContainer;

extern Timer GlobalTimer;

void
CudaLoop::BasicInit(const ParamFile& paramfile,
                    const ProblemContainer* PC,
                    const FunctionalContainer* FC)
{
  GetMultiLevelSolverPointer() = new CudaMultiLevelSolver;
  StdLoop::BasicInit(paramfile, PC, FC);
}

} // namespace Gascoigne
