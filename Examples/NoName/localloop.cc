#include  "localloop.h"
#include  "localmeshagent.h"

#include  "compose_name.h"
#include  "backup.h"
#include  "adaptordata.h"
#include  "diplomantenadaptor.h"
#include  "malteadaptor.h"
#include  "filescanner.h"
#include  "monitoring.h"
#include  <iomanip>
#include  "gostream.h"
#include  "localmultilevelsolver.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

void LocalLoop::BasicInit(const ParamFile* paramfile)
{
  GetMeshAgentPointer() = new LocalMeshAgent;
  GetMeshAgent()->BasicInit(paramfile);
  GetMultiLevelSolverPointer() = new LocalMultiLevelSolver;

  StdLoop::BasicInit(paramfile);
}
