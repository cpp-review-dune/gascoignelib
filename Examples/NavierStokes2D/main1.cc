#include  "problemdescriptor1.h"
#include  "stdloop.h"
#include  "starter.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  Starter S(argc, argv, "bench.param");

  ProblemDescriptor1 LPD;
  LPD.BasicInit(S.GetParamFile());

  StdLoop loop;
  loop.BasicInit(S.GetParamFile());
  loop.run(&LPD);

  return 0;
}
