#include  "problemdescriptorinitial.h"
#include  "localequation.h"
#include  "localinitialcondition.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptorInitial::ConstructEquation()
{
  GetEquationPointer() = new LocalEquation(GetParamFile());
}

/* ----------------------------------------- */

void ProblemDescriptorInitial::ConstructInitialCondition()
{
  const LocalEquation* LEQ = dynamic_cast<const LocalEquation*>(GetEquation());
  GetInitialConditionPointer() = new LocalInitialCondition(LEQ);
}
