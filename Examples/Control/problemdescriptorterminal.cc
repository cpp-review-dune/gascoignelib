#include  "problemdescriptorterminal.h"
#include  "localequation.h"
#include  "localterminalcondition.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptorTerminal::ConstructEquation()
{
  GetEquationPointer() = new LocalEquation(GetParamFile());
}

/* ----------------------------------------- */

void ProblemDescriptorTerminal::ConstructInitialCondition()
{
  GetInitialConditionPointer() = new LocalTerminalCondition();
}
