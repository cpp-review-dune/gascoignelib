#include  "problemdescriptor2.h"
#include  "stdloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh2.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor2 LPD;
  LPD.BasicInit(&paramfile);

  StdLoop loop;
  loop.BasicInit(&paramfile);
  loop.run(&LPD);

  return 0;
}
