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
  loop.init("Results/forward",0,&LPD);

  return 0;
}
