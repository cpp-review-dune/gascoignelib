#include  "local1.h"
#include  "stdloop.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("box.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

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
