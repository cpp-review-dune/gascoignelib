#include  "problemdescriptor1.h"
#include  "laplace2d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "polynomialexactsolution.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructEquation()
{
  GetEquationPointer() = new Laplace2d;
}

/* ----------------------------------------- */

void ProblemDescriptor1::ConstructExactSolution()
{
  GetExactSolutionPointer() = new PolynomialExactSolution();
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
