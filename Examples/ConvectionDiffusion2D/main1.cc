#include  "localloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh1.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  LocalLoop loop;
  loop.BasicInit(&paramfile);
  loop.run();

  return 0;
}
