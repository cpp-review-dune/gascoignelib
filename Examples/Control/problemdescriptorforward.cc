#include  "problemdescriptorforward.h"
#include  "localequation.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptorForward::ConstructEquation()
{
  GetEquationPointer() = new LocalEquation(GetParamFile());
}
