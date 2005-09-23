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

  virtual DiscretizationInterface* CreateOtherDiscretization() const;

  double ScalarProduct(nvector<double>& eta, const GlobalVector& f, 
		       const GlobalVector& z) const;

  double ScalarProduct(nvector<double>& eta, const VectorInterface& gf, 
		       const VectorInterface& gz) const;
  double ScalarProductWithFluctuations(nvector<double>& eta, 
				       const VectorInterface& gf, 
				       const VectorInterface& gz) const;

  void PrimalResidualsHigher(VectorInterface& gf, 
			     const VectorInterface& gu);

  void DualResidualsHigher(VectorInterface& gf, const VectorInterface& gu, 
			   const VectorInterface& gz, 
			   const ProblemDescriptorInterface& PDI);

 public:

  DwrQ1Q2(SolverInterface& SR);
  virtual ~DwrQ1Q2() {};

  double Estimator(nvector<double>& eta, VectorInterface& gf, 
		   const VectorInterface& gu, const VectorInterface& gz,
		   const ProblemDescriptorInterface& PDI);
};
}
/*-------------------------------------------------------*/

#endif
