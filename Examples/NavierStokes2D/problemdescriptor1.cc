#include  "problemdescriptor1.h"
#include  "navierstokesgls.h"
#include  "benchmarkdirichletdata.h"
#include  "zerorighthandsidedata.h"
#include  "zeroexactsolution.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructEquation()
{
  GetEquationPointer() = new NavierStokesGls(GetParamFile());
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructDirichletData()
{
  GetDirichletDataPointer() = new BenchMarkDirichletData(GetParamFile());
}
