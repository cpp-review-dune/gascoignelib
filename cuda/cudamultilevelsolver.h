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
public:
  CudaMultiLevelSolver();
  CudaMultiLevelSolver(const MeshAgent* GMGM,
                       const ParamFile& paramfile,
                       const ProblemContainer* PC,
                       const FunctionalContainer* FC = NULL);
  virtual ~CudaMultiLevelSolver();

  StdSolver* NewSolver(IndexType) override;
  void NewMgInterpolator() override;

  void LinearMg(IndexType finelevel,
                IndexType coarselevel,
                const Matrix& A,
                Vector& u,
                const Vector& f,
                CGInfo& info) override;
  void mgstep(std::vector<double>& res,
              std::vector<double>& rw,
              IndexType l,
              IndexType maxl,
              IndexType minl,
              std::string& p0,
              std::string p,
              const Matrix& A,
              Vector& u,
              Vector& b,
              Vector& v) override;
};

} // namespace Gascoigne

/*----------------------------   solver.h     ---------------------------*/
/* end of #ifndef cudamultilevelsolver_H */
#endif
/*----------------------------   solver.h     ---------------------------*/
