#include "solvers.h"


using namespace std;

extern double __DT;


namespace Gascoigne
{

  double NSSolver::ComputeResidualFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, 
					     const ResidualFunctional* FP) const 
  {
    double j = StdSolver::ComputeResidualFunctional(gf,gu,gz,FP);
    return j/__DT;
  }

}
