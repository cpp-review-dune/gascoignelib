#include  "local.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("bench.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor1 LPD;
  LPD.BasicInit(&paramfile);

  LocalLoop loop;
  loop.BasicInit(&paramfile);
  loop.run(&LPD);

  return 0;
}
