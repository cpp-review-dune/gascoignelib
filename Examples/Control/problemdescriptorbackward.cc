#include  "problemdescriptorbackward.h"
#include  "backwardequation.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptorBackward::ConstructEquation()
{
  GetEquationPointer() = new BackwardEquation(GetParamFile());
}
