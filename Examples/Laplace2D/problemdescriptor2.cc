#include  "problemdescriptor2.h"
#include  "laplace.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "polynomialexactsolution.h"
#include  "neumanndatabyexactsolution.h"

using namespace std;

/* ----------------------------------------- */

void ProblemDescriptor2::ConstructEquation()
{
  GetEquationPointer() = new Laplace;
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
  cerr << "hallooo      ProblemDescriptor2::ConstructNeumannData()      \n";
  GetNeumannDataPointer() = new NeumannDataByExactSolution(GetEquation(),GetExactSolution());
}
