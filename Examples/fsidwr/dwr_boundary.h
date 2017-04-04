#ifndef __DWR_BOUNDARY_H__
#define __DWR_BOUNDARY_H__

#include "solverinterface.h"
#include "stdmultilevelsolver.h"
#include "dwr_discretization.h"
#include "alediscretization.h"
#include <set>
#include <map>


/*-------------------------------------------------------*/

namespace Gascoigne
{
  class DWRBoundary
  {
  protected:

    const ParamFile* _paramfile;
    
    
    StdMultiLevelSolver*       __MS;
    const std::string&         __primalproblem,__dualproblem;
    
    SolverInterface*           __S;
    DiscretizationInterface*   __D;

    DiscretizationInterface    *__higher, *__lower;
    
    double ScalarProduct(nvector<double>& eta, const GlobalVector& f, const GlobalVector& z) const;

    double ScalarProduct(nvector<double>& eta, const VectorInterface& gf,  const VectorInterface& gz) const;
    
  public:

    DWRBoundary(MultiLevelSolverInterface* SR, const std::string primalproble, const std::string dualproblem, const ParamFile* __pf);
      virtual ~DWRBoundary() {};
    
    double ScalarProductWithFluctuations(nvector<double>& eta, 
					 const VectorInterface& gf, 
					 const VectorInterface& gz) const;
    
    void ResidualsHigher(VectorInterface& gf, const VectorInterface& gu, DiscretizationInterface* D, double s);
    
    double Estimator(nvector<double>& eta, VectorInterface& gf, 
		     VectorInterface& gu, VectorInterface& gz);
    double EstimatorEasy(nvector<double>& eta, VectorInterface& gf, 
			 VectorInterface& gu, VectorInterface& gz);
    double EstimatorAdjoint(nvector<double>& eta, VectorInterface& gf, 
			    VectorInterface& gu, VectorInterface& gz);



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
