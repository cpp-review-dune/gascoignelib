#include  "problemdescriptorinitial.h"
#include  "localtimeloop.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("initial.param");

  ProblemDescriptorInitial LPD;
  LPD.BasicInit(&paramfile);

  LocalTimeLoop loop;
  loop.BasicInit(&paramfile);
  loop.NewMesh(&LPD);
  loop.init("Results/forward",0);

  return 0;
}
