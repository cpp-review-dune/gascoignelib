#include "solvers.h"

extern double CHIGAUSS;
extern bool   PRIMALPROBLEM;

namespace Gascoigne{

void MySolver::AssembleMatrix(const VectorInterface& gu, double d)
{
  if (PRIMALPROBLEM)
    {
      CHIGAUSS = 0.5-0.5*sqrt(1.0/3.0);
      StdSolver::AssembleMatrix(gu,d*0.5);
      CHIGAUSS = 0.5+0.5*sqrt(1.0/3.0);
      StdSolver::AssembleMatrix(gu,d*0.5);
    }
  else
    StdSolver::AssembleMatrix(gu,d);
}

void MySolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  if (PRIMALPROBLEM)
    {
      CHIGAUSS = 0.5-0.5*sqrt(1.0/3.0);
      StdSolver::Form(gy,gx,d*0.5);
      CHIGAUSS = 0.5+0.5*sqrt(1.0/3.0);
      StdSolver::Form(gy,gx,d*0.5);
    }
  else
    StdSolver::Form(gy,gx,d);
}
  
}
