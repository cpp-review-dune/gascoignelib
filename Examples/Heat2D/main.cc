#include  "local.h"
#include  "stdtimeloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("heat", &LPD);
  
  StdTimeLoop loop;

  loop.BasicInit(&paramfile,&PC);
  loop.run("heat");

  return 0;
}
