#ifndef cudamultilevelsolver_H
#define cudamultilevelsolver_H

#include <string>
#include <vector>

#include "stdmultilevelsolver.h"
#include "stdsolver.h"
#include "vectorinterface.h"

namespace Gascoigne {
class CGInfo;
class StdSolver;

class CudaMultiLevelSolver : public virtual StdMultiLevelSolver
{
protected:
  void ActivateCuda(std::initializer_list<const Vector*> vectors);
  void DeactivateCuda(std::initializer_list<Vector*> vectors);

public:
  CudaMultiLevelSolver();
  CudaMultiLevelSolver(const MeshAgent* GMGM,
                       const ParamFile& paramfile,
                       const ProblemContainer* PC,
                       const FunctionalContainer* FC = NULL);
  virtual ~CudaMultiLevelSolver();

  StdSolver* NewSolver(int) override;
  void NewMgInterpolator() override;

  void NewtonLinearSolve(const Matrix& A,
                         Vector& x,
                         const Vector& b,
                         CGInfo& info);

  void RestrictZero(IndexType level, Vector& b, const Vector& v) const override;
  void ProlongateAdd(IndexType level,
                     Vector& b,
                     const Vector& v) const override;
};

} // namespace Gascoigne

/*----------------------------   solver.h     ---------------------------*/
/* end of #ifndef cudamultilevelsolver_H */
#endif
/*----------------------------   solver.h     ---------------------------*/
