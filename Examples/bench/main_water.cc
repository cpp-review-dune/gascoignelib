#include  "problem.h"
#include  "loop.h"

using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  // parameter file as argument
  ParamFile paramfile("gascoigne.param");
  if(argc>=2) 
    paramfile.SetName(argv[1]);

  WaterProblem WPD;
  WPD.BasicInit(&paramfile);


  ProblemContainer PC;
  PC.AddProblem("water",  &WPD);

  // output functionals (not in use now)
  FunctionalContainer FC;

  // loop for program control
  Loop loop;
  loop.BasicInit(&paramfile, &PC, &FC);

  // start solving
  loop.runwater("water");

  return 0;
}
