#ifndef __EnergyEstimator_h

#include "solverinterface.h"
#include "meshinterpretorinterface.h"

/*-------------------------------------------------------*/

namespace Gascoigne
{

class EnergyEstimator
{
 protected:

  SolverInterface& S;
  const ProblemDescriptorInterface* primalproblem;
  MeshInterpretorInterface*         discretization;

 public:

  EnergyEstimator(SolverInterface& SR);
  ~EnergyEstimator() {};

  double Estimator(nvector<double>& eta, BasicGhostVector& gu, 
		   const BasicGhostVector& gf);
};

/*-------------------------------------------------------*/

}
#endif
