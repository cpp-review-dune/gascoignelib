#include  "problemdescriptor.h"
#include  "localequation.h"
#include  "localdirichletdata.h"
#include  "zerorighthandsidedata.h"
#include  "zeroexactsolution.h"

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
