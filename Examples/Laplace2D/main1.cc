#include  "local1.h"
#include  "stdloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh1.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  /////////////
  // Loop
  /////////////
  StdLoop loop;
  loop.BasicInit(&paramfile);

  /////////////
  // Functionals
  /////////////
  LocalDragFunctional   j0; 
  LocalDomainFunctional j1;
  loop.AddFunctional(&j0);
  loop.AddFunctional(&j1);
  
  loop.run(&LPD);

  return 0;
}
