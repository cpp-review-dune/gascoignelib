#include  "problemdescriptor.h"
#include  "localloop.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  std::string paramfile("bench.param");
  if(argc>=2) {
    paramfile = argv[1];
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(paramfile);

  LocalLoop loop;
  loop.BasicInit(paramfile,&LPD);
  loop.run();

  return 0;
}
