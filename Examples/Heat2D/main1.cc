#include  "problemdescriptor1.h"
#include  "stdtimeloop.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  std::string paramfile("mesh.param");
  if(argc>=2) {
    paramfile = argv[1];
  }

  ProblemDescriptor1 LPD;
  LPD.BasicInit(paramfile);

  StdTimeLoop loop;
  loop.BasicInit(paramfile);
  loop.run(&LPD);

  return 0;
}
