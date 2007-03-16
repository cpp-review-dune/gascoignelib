#include  "local.h"


using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC) 
  {
    GetMeshAgentPointer() = new BenchMarkMeshAgent;
    StdLoop::BasicInit(paramfile, PC, FC);
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

  // Functionals
  FunctionalContainer FC;
  LocalDomainFunctionals_FlagForce df_flagforce_lift("lift");
  LocalDomainFunctionals_FlagForce df_flagforce_drag("drag");
  FC.AddFunctional("lift",&df_flagforce_lift);
  FC.AddFunctional("drag",&df_flagforce_drag);
  
  LocalLoop loop;
  loop.BasicInit(&paramfile, &PC, &FC);
  loop.run("navier stokes");

  return 0;
}
