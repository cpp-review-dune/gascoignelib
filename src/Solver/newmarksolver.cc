#include "newmarksolver.h"

using namespace Gascoigne;
using namespace std;

/*----------------------------------------------------------------------------*/

void NewmarkSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  FormWithoutMass(gy,gx,d);
  MassMatrixVector(gy,gx,d);
}  

/*----------------------------------------------------------------------------*/

void NewmarkSolver::FormWithoutMass(VectorInterface& gy, const VectorInterface& gx, double d, double s) const
{
  HNAverage(gx);
  HNAverageData();
  
  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  double alpha = 0.25*dt*dt;
  
  GetDiscretization()->Form(GetGV(gy),GetGV(gx),*EQ,d*alpha);
  
  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
    {
      const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
      GetDiscretization()->BoundaryForm(GetGV(gy),GetGV(gx),BM->GetBoundaryEquationColors(),*BE,d*sqrt(alpha)*s);
    }
  HNZero(gx);
  HNZeroData();
  HNDistribute(gy);
  SubtractMeanAlgebraic(gy);
}

/*----------------------------------------------------------------------------*/

void NewmarkSolver::AssembleMatrix(const VectorInterface& gu, double d)
{
  assert(GetMatrix());
  
  double alpha = 0.25*dt*dt;
  
  const GlobalVector& u = GetGV(gu);
  HNAverage(gu);
  HNAverageData();
  
  GetDiscretization()->Matrix(*GetMatrix(),u,*GetProblemDescriptor()->GetEquation(),d*alpha);
  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
    {
      const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
      GetDiscretization()->BoundaryMatrix(*GetMatrix(),u,BM->GetBoundaryEquationColors(),*BE,d*sqrt(alpha));
    }
  DirichletMatrix();
  HNZero(gu);
  HNZeroData();
  GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),_TP,d);
  StdSolver::DirichletMatrix();
}

/*----------------------------------------------------------------------------*/

