#include  "local.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile) 
    {
      GetMeshAgentPointer() = new BenchMarkMeshAgent;
      StdLoop::BasicInit(paramfile);
    }
};

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("bench.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  LocalLoop loop;
  loop.BasicInit(&paramfile);
  loop.run(&LPD);

  return 0;
}
