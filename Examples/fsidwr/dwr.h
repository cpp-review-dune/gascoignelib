#ifndef __DWR_h

#include "solverinterface.h"
#include "alediscretization.h"
/*-------------------------------------------------------*/

namespace Gascoigne
{
  class Dwr
  {
  protected:

    SolverInterface& S;
    const ProblemDescriptorInterface* primalproblem;
    const ProblemDescriptorInterface* dualproblem;
    DiscretizationInterface*         discretization;
    DiscretizationInterface*         other;

    virtual DiscretizationInterface* CreateOtherDiscretization() const;

    double ScalarProduct(nvector<double>& eta, const GlobalVector& f, 
			 const GlobalVector& z) const;

    double ScalarProduct(nvector<double>& eta, const VectorInterface& gf, 
			 const VectorInterface& gz) const;

  public:

    Dwr(SolverInterface& SR, 
	const ProblemDescriptorInterface* primal,
	const ProblemDescriptorInterface* dual);

    double ScalarProductWithFluctuations(nvector<double>& eta, 
					 const VectorInterface& gf, 
					 const VectorInterface& gz) const;
    double ScalarProductWithFluctuationsPrimal(nvector<double>& eta, 
					       const VectorInterface& gf, 
					       const VectorInterface& gz) const;

    void PrimalResidualsHigher(VectorInterface& gf, 
			       const VectorInterface& gu);

    void DualResidualsHigher(VectorInterface& gf, const VectorInterface& gu, 
			     const VectorInterface& gz, 
			     const ProblemDescriptorInterface& PDI);

    double Estimator(nvector<double>& eta, VectorInterface& gf, 
		     VectorInterface& gu, VectorInterface& gz);
    double EstimatorEnergy(nvector<double>& eta, VectorInterface& gf,
			   const VectorInterface& gu);

    void ReInitInterface(DiscretizationInterface* DISC);
    void reinit_element_2d(int en, const nvector<int>& indices,
			   HASHMAP<int, std::vector<int> >& solid_interface_cells, 
			   HASHMAP<int, std::vector<int> >& fluid_interface_cells,
			   HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
			   HASHSET<int> & interface_nodes,
			   std::set<int>& fluid_nodes, std::set<int>& solid_nodes);
    void reinit_element_3d(int en, const nvector<int>& indices, 
			   HASHMAP<int, std::vector<int> >& solid_interface_cells, 
			   HASHMAP<int, std::vector<int> >& fluid_interface_cells,
			   HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
			   HASHSET<int> & interface_nodes,
			   std::set<int>& fluid_nodes, std::set<int>& solid_nodes);

  };
}
/*-------------------------------------------------------*/

#endif
