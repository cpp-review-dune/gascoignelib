#include  "problemdescriptor1.h"
#include  "stdtimeloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor1 LPD;
  LPD.BasicInit(&paramfile);

  StdTimeLoop loop;

  loop.BasicInit(&paramfile);
  loop.run(&LPD);

  return 0;
}
