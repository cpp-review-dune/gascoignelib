#include  "problemdescriptor1.h"
#include  "localequation.h"
#include  "localinitialcondition.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructEquation()
{
  GetEquationPointer() = new LocalEquation(GetParamFile());
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructInitialCondition()
{
  const LocalEquation* LEQ = dynamic_cast<const LocalEquation*>(GetEquation());
  GetInitialConditionPointer() = new LocalInitialCondition(LEQ);
}
