#include  "problemdescriptor1.h"
#include  "zerorighthandsidedata.h"
#include  "localequation.h"
#include  "localdirichletdata.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructEquation()
{
  GetEquationPointer() = new LocalEquation(GetParamFile());
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructRightHandSideData()
{
  int ncomp = 1;
  GetRightHandSideDataPointer() = new ZeroRightHandSideData(ncomp);
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructDirichletData()
{
  GetDirichletDataPointer() = new LocalDirichletData;
}
