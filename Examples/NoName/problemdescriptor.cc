#include  "problemdescriptor.h"
#include  "localequation.h"
#include  "localdirichletdata.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor::ConstructEquation()
{
  GetEquationPointer() = new LocalEquation(GetParamFile());
}

/* ----------------------------------------- */

void ProblemDescriptor::ConstructDirichletData()
{
  GetDirichletDataPointer() = new LocalDirichletData(GetParamFile());
}
