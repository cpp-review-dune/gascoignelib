#include  "problemdescriptorinitial.h"
#include  "localtimeloop.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  std::string paramfile("initial.param");
  if(argc>=2) {
    paramfile = argv[1];
  }

  ProblemDescriptorInitial LPD;
  LPD.BasicInit(paramfile);

  LocalTimeLoop loop;
  loop.BasicInit(paramfile,&LPD);
  loop.NewMesh();
  loop.init("Results/forward",0);

  return 0;
}
