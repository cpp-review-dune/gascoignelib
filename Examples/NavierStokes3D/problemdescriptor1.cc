#include  "problemdescriptor1.h"
#include  "stokesgls3d.h"
#include  "navierstokesgls3d.h"
#include  "benchmark3ddirichletdata.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructEquation()
{
  GetEquationPointer() = new NavierStokesGls3d(GetParamFile());
//   GetEquationPointer() = new StokesGls3d(_paramfile);
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructDirichletData()
{
  GetDirichletDataPointer() = new BenchMark3dDirichletData();
}
