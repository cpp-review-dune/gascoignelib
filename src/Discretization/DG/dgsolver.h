
/*----------------------------   dgsolver.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dgsolver_H
#define __dgsolver_H
/*----------------------------   dgsolver.h     ---------------------------*/

#include "stdsolver.h"
#include "stdmultilevelsolver.h"

namespace Gascoigne
{
  

class DGSolver : public StdSolver
{
public:

  std::string GetName() const
  {
    return "DG Solver";
  }
  
  void
  Visu(const std::string &name, const VectorInterface &gu, int i) const;
};

  class DGMultiLevelSolver : public StdMultiLevelSolver
  {
  public:
    SolverInterface* NewSolver(int solverlevel)
    {
      return new DGSolver;
    }
  };
  
}



/*----------------------------   dgsolver.h     ---------------------------*/
/* end of #ifndef __dgsolver_H */
#endif
/*----------------------------   dgsolver.h     ---------------------------*/
