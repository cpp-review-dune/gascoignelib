#include  "localstokes.h"
#include  "localloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{

  ParamFile paramfile("periodic.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  LocalStokesProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("stokes", &LPD);
    
  FunctionalContainer FC;
  LocalDragFunctional LDF;
  LocalLiftFunctional LLF;
  FC.AddFunctional("Drag", &LDF);
  FC.AddFunctional("Lift", &LLF);

  LocalLoop loop;

  loop.BasicInit(&paramfile,&PC,&FC);
  loop.run("stokes");

  return 0;
}
