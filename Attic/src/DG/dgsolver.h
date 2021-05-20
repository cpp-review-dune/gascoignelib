
/*----------------------------   dgsolver.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dgsolver_H
#define __dgsolver_H
/*----------------------------   dgsolver.h     ---------------------------*/

#include "stdmultilevelsolver.h"
#include "stdsolver.h"

namespace Gascoigne {

class DGSolver : public StdSolver
{
public:
  std::string GetName() const { return "DG Solver"; }

  void Visu(const std::string& name, const Vector& gu, int i) const;
};

class DGMultiLevelSolver : public StdMultiLevelSolver
{
public:
  StdSolver* NewSolver(int solverlevel) { return new DGSolver; }

  // Mehrgitter-Transfer nicht implementiert, alles nur auf dem feisnten gitter

  void Transfer(int high, int low, Vector& u) const { assert(0); }
  void AddNodeVector(const std::string& name, Vector& gq)
  {
    GetSolver()->AddNodeVector(name, gq);
  }
  void DeleteNodeVector(const std::string& name)
  {
    GetSolver()->DeleteNodeVector(name);
  }

  void AssembleMatrix(Vector& u)
  {
    GetSolver()->MatrixZero();
    GetSolver()->AssembleMatrix(u, 1.);
  }
  void ComputeIlu(Vector& u) { GetSolver()->ComputeIlu(u); }
};

} // namespace Gascoigne

/*----------------------------   dgsolver.h     ---------------------------*/
/* end of #ifndef __dgsolver_H */
#endif
/*----------------------------   dgsolver.h     ---------------------------*/
