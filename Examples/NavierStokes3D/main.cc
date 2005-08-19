#include  "local.h"
#include  "stdloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("benchcircle.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("navierstokes", &LPD);

  StdLoop loop;

  loop.BasicInit(&paramfile,&PC);
  loop.run("navierstokes");

  return 0;
}
