#include  "problemdescriptor.h"
#include  "localloop.h"
#include  "starter.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  Starter S(argc, argv, "bench.param");

  ProblemDescriptor LPD;
  LPD.BasicInit(S.GetParamFile());

  LocalLoop loop;
  loop.BasicInit(S.GetParamFile(),LPD);
  loop.run();

  return 0;
}
