#include  "problemdescriptor1.h"
#include  "laplace3d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "polynomialexactsolution3d.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructEquation()
{
  GetEquationPointer() = new Laplace3d(GetParamFile());
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructExactSolution()
{
  GetExactSolutionPointer() = new PolynomialExactSolution3d();
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructRightHandSideData()
{
  GetRightHandSideDataPointer() = new RightHandSideDataByEquation(GetEquation(), GetExactSolution());
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructDirichletData()
{
  GetDirichletDataPointer() = new DirichletDataByExactSolution(GetExactSolution());
}
