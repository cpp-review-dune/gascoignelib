#include  "localloop.h"
#include  "starter.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  Starter S(argc, argv, "mesh1.param");

  LocalLoop loop;
  loop.BasicInit(S.GetParamFile());
  loop.run();

  return 0;
}
