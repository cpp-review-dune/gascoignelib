#include  "problemdescriptor2.h"
#include  "laplace2d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "polynomialexactsolution.h"
#include  "neumanndatabyexactsolution.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor2::ConstructEquation()
{
  GetEquationPointer() = new Laplace2d;
}

/* ----------------------------------------- */

void ProblemDescriptor2::ConstructExactSolution()
{
  GetExactSolutionPointer() = new PolynomialExactSolution();
}

/* ----------------------------------------- */

void ProblemDescriptor2::ConstructRightHandSideData()
{
  GetRightHandSideDataPointer() = new RightHandSideDataByEquation(GetEquation(), GetExactSolution());
}

/* ----------------------------------------- */

void ProblemDescriptor2::ConstructDirichletData()
{
  GetDirichletDataPointer() = new DirichletDataByExactSolution(GetExactSolution());
}

/* ----------------------------------------- */

void ProblemDescriptor2::ConstructNeumannData()
{
  GetNeumannDataPointer() = new NeumannDataByExactSolution(GetEquation(),GetExactSolution());
}
