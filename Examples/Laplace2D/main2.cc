#include  "problemdescriptor2.h"
#include  "stdloop.h"
#include  "starter.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  Starter S(argc, argv, "mesh2.param");

  ProblemDescriptor2 LPD;
  LPD.BasicInit(S.GetParamFile());

  StdLoop loop;
  loop.BasicInit(S.GetParamFile());
  loop.run(&LPD);

  return 0;
}
