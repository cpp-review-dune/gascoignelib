/*----------------------------   solvers.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solvers_H
#define __solvers_H
/*----------------------------   solvers.h     ---------------------------*/


#include "stdmultilevelsolver.h"
#include "stdsolver.h"


using namespace std;

namespace Gascoigne
{

  class NSSolver : public StdSolver
  {
  private:
  public:
    std::string GetName() const {return "FSI Solver";}
    double ComputeResidualFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const ResidualFunctional* FP) const;
    
    
  };
  
  
  
  class MLS : public StdMultiLevelSolver
  {
  public:
    std::string GetName() const {return "NS MultiLevelSolver";}
    
    SolverInterface* NewSolver(int solverlevel)
    { return new NSSolver; }
    
    const NSSolver* GetNSSolver(int l) const
    {
      assert(dynamic_cast<const NSSolver* > (GetSolver(l)));
      return dynamic_cast<const NSSolver* > (GetSolver(l));
    }
  };
  
}

/*----------------------------   solvers.h     ---------------------------*/
/* end of #ifndef __solvers_H */
#endif
/*----------------------------   solvers.h     ---------------------------*/
