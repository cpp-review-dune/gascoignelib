#ifndef __DwrQ1Q2_h

#include "solverinterface.h"

/*-------------------------------------------------------*/

namespace Gascoigne
{
class DwrQ1Q2
{
 protected:

  SolverInterface& S;
  const ProblemDescriptorInterface* primalproblem;
  DiscretizationInterface*         discretization;

  double ScalarProduct(nvector<double>& eta, const GlobalVector& f, 
		       const GlobalVector& z) const;

  double ScalarProduct(nvector<double>& eta, const BasicGhostVector& gf, 
		       const BasicGhostVector& gz) const;
  double ScalarProductWithFluctuations(nvector<double>& eta, 
				       const BasicGhostVector& gf, 
				       const BasicGhostVector& gz) const;

  void PrimalResidualsHigher(BasicGhostVector& gf, 
			     const BasicGhostVector& gu);

  void DualResidualsHigher(BasicGhostVector& gf, const BasicGhostVector& gu, 
			   const BasicGhostVector& gz, 
			   const ProblemDescriptorInterface& PDI);

 public:

  DwrQ1Q2(SolverInterface& SR);
  ~DwrQ1Q2() {};

  double Estimator(nvector<double>& eta, BasicGhostVector& gf, 
		   const BasicGhostVector& gu, const BasicGhostVector& gz,
		   const ProblemDescriptorInterface& PDI);
};
}
/*-------------------------------------------------------*/

#endif
