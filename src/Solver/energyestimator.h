#ifndef __EnergyEstimator_h

#include "solverinterface.h"
#include "q1.h"

/*-------------------------------------------------------*/

namespace Gascoigne
{

class EnergyEstimator
{
 protected:

  SolverInterface& S;
  const ProblemDescriptorInterface* primalproblem;
  Q1*                               discretization;
  double                            _d_visc;
  std::string                       _s_energytype;

 public:

  EnergyEstimator(SolverInterface& SR);
  ~EnergyEstimator() {};

  double Estimator(DoubleVector& eta, VectorInterface& gu, 
		   const VectorInterface& gf);
};

/*-------------------------------------------------------*/

}
#endif
