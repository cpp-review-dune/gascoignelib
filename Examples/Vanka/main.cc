
#include "local.h"
#include "loop.h"

using namespace Gascoigne;
using namespace std;

int
main(int argc, char** argv)
{
  ParamFile pf("fsi-3.param");
  if (argc == 2)
    pf.SetName(argv[1]);

  ProblemDescriptor<2> Problem2d;
  Problem2d.BasicInit(&pf);

  ProblemContainer PC2d;
  PC2d.AddProblem("fsi", &Problem2d);

  FunctionalContainer FC2d;
  Loop<2> loop;

  loop.BasicInit(&pf, &PC2d, &FC2d);
  loop.run("fsi");

  return 0;
}
