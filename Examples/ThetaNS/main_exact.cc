
#include  "local.h"
#include  "loop.h"
#include "functionals.h"


using namespace Gascoigne;
using namespace std;




/*---------------------------------------------------*/



int main(int argc, char** argv)
{
  ParamFile pf("bench.param");
  assert(argc==1);

  // FUNCTIONAL

  FunctionalContainer FC2d;
  VMEAN vmean;
  FC2d.AddFunctional("0 vmean", &vmean);

  Blift blift(&pf);
  Bdrag bdrag(&pf);
  Pdrop pdrop;
  FC2d.AddFunctional("1 blift", &blift);
  //  FC2d.AddFunctional("2 bdrag", &bdrag);
  //  FC2d.AddFunctional("4 pdrop", &pdrop);

  // PROBLEM
  
  ProblemDescriptor2d primal;
  primal.BasicInit(&pf);

  ProblemContainer PC2d;
  PC2d.AddProblem("primal", &primal);

  Loop loop;

  loop.BasicInit(&pf,&PC2d,&FC2d);
  loop.run_exact_fst("primal");

  return 0;
}
