#include  "local.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  /////////////
  // Equation
  /////////////
//  ProblemDescriptor LPD;
//  LPD.BasicInit(&paramfile);
  ProjectionProblemDescriptor PD;
  PD.BasicInit(&paramfile);

  /////////////
  // Loop
  /////////////
  LocalLoop loop;
  loop.BasicInit(&paramfile);

  loop.run(&PD);

  return 0;
}
