#include  "problemdescriptorinitial.h"
#include  "localtimeloop.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("initial.param");

  ProblemDescriptorInitial LPD;
  LPD.BasicInit(&paramfile);

  LocalTimeLoop loop;
  loop.BasicInit(&paramfile,&LPD);
  loop.NewMesh();
  loop.init("Results/forward",0);

  return 0;
}
