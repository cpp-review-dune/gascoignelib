#ifndef __EnergyEstimator_h

#include "solverinterface.h"

/*-------------------------------------------------------*/

namespace Gascoigne
{

class EnergyEstimator
{
 protected:

  SolverInterface& S;
  const ProblemDescriptorInterface* primalproblem;
  DiscretizationInterface*          discretization;

 public:

  EnergyEstimator(SolverInterface& SR);
  ~EnergyEstimator() {};

  double Estimator(nvector<double>& eta, BasicGhostVector& gu, 
		   const BasicGhostVector& gf);
};

/*-------------------------------------------------------*/

}
#endif
