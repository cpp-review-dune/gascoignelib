#include  "problemdescriptor1.h"
#include  "stdloop.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  std::string paramfile("box.param");
  if(argc>=2) {
    paramfile = argv[1];
  }

  ProblemDescriptor1 LPD;
  LPD.BasicInit(paramfile);

  StdLoop loop;
  loop.BasicInit(paramfile);
  loop.run(&LPD);

  return 0;
}
