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

/* ----------------------------------------- */

void LocalLoop::BasicInit(const string& pfile, const ProblemDescriptorInterface* PD)
{
  GetMeshAgentPointer() = new LocalMeshAgent;
  GetMultiLevelSolverPointer() = new LocalMultiLevelSolver;

  StdLoop::BasicInit(pfile,PD);
}
