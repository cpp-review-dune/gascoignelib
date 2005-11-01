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

  ProblemContainer PC;
  PC.AddProblem("laplace3d", &LPD);

  FunctionalContainer FC;
  LocalDragFunctional   j0; 
  LocalDomainFunctional j1;
  FC.AddFunctional("drag", &j0);
  FC.AddFunctional("domain", &j1);
  
  
  StdLoop loop;

  loop.BasicInit(&paramfile, &PC, &FC);

  loop.run("laplace3d");

  return 0;
}
