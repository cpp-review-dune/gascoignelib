#include  "local.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC) 
    {
      GetMeshAgentPointer() = new BenchMarkMeshAgent;
      StdLoop::BasicInit(paramfile, PC);
    }
};

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("bench.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("navier stokes", &LPD);
  
  LocalDomainFunctionals_FlagForce df_flagforce_lift("lift");
  LocalDomainFunctionals_FlagForce df_flagforce_drag("drag");

  LocalLoop loop;
  loop.BasicInit(&paramfile, &PC);
  loop.AddFunctional(&df_flagforce_lift);
  loop.AddFunctional(&df_flagforce_drag);
  loop.run("navier stokes");

  return 0;
}
