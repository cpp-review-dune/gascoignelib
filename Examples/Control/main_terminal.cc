#include  "problemdescriptorterminal.h"
#include  "localtimeloop.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("initial.param");
  if(argc!=3) {
    cerr << "usage: \"terminal: name last\"";
  }
  string  name(argv[1]);
  int last = atoi(argv[2]);


  ProblemDescriptorTerminal LPD;
  LPD.BasicInit(&paramfile);
  LPD.GetInitialCondition();

  LocalTimeLoop loop;
  loop.BasicInit(&paramfile);
  loop.ReInit(&LPD);
  loop.AddNodeVector(name);
  loop.init("Results/backward",last,&LPD);

  return 0;
}
