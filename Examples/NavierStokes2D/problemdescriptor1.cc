#include  "problemdescriptor1.h"
#include  "navierstokesgls2d.h"
#include  "benchmarkdirichletdata.h"
#include  "zerorighthandsidedata.h"
#include  "zeroexactsolution.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructEquation()
{
  GetEquationPointer() = new NavierStokesGls2d(GetParamFile());
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructDirichletData()
{
  GetDirichletDataPointer() = new BenchMarkDirichletData();
}
