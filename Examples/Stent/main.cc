#include  "local.h"
#include  "loop.h"

using namespace Gascoigne;
using namespace std;


int main(int argc, char** argv)
{
  ParamFile pf("stent.param");
  if (argc==2)
    pf.SetName(argv[1]);

  
  PD Problem;
  Problem.BasicInit(&pf);

  ProblemContainer PC2d;
  PC2d.AddProblem("fsi", &Problem);

  FunctionalContainer FC2d;

  Loop loop;

  loop.BasicInit(&pf,&PC2d,&FC2d);
  loop.run("fsi");

  return 0;
}
