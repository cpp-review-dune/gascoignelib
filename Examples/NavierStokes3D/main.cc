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

  StdLoop loop;

  loop.BasicInit(&paramfile);
  loop.run(&LPD);

  return 0;
}