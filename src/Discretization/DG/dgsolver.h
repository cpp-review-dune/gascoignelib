
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
    StdSolver* NewSolver(int solverlevel)
    {
      return new DGSolver;
    }
    
    
    
    // Mehrgitter-Transfer nicht implementiert, alles nur auf dem feisnten gitter
    
    
    void Transfer(int high, int low, VectorInterface& u) const
    { assert(0); }
    void AddNodeVector(const std::string& name, VectorInterface& gq)
    {
      GetSolver()->AddNodeVector(name,gq);
    }
    void DeleteNodeVector(const std::string& name)
    {
      GetSolver()->DeleteNodeVector(name);
    }
    
    void AssembleMatrix(VectorInterface& u)
    {
      GetSolver()->MatrixZero();
      GetSolver()->AssembleMatrix(u,1.);
    }
    void ComputeIlu(VectorInterface& u)
    {
      GetSolver()->ComputeIlu(u);
    }
  };
  
}



/*----------------------------   dgsolver.h     ---------------------------*/
/* end of #ifndef __dgsolver_H */
#endif
/*----------------------------   dgsolver.h     ---------------------------*/
