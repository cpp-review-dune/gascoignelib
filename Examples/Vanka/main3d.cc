#include "local.h"
#include "loop.h"

using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

int
main(int argc, char** argv)
{
  ParamFile pf("box3d.param");
  if (argc == 2)
    pf.SetName(argv[1]);

  ProblemDescriptor<3> Problem3d;
  Problem3d.BasicInit(&pf);

  ProblemContainer PC3d;
  PC3d.AddProblem("fsi", &Problem3d);

  FunctionalContainer FC3d;

  Loop<3> loop;

  loop.BasicInit(&pf, &PC3d, &FC3d);
  loop.run("fsi");

  return 0;
}
