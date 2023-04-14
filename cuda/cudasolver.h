/*----------------------------   cudasolver.h     ---------------------------*/
/*      $Id:$                 */
#ifndef cudasolver_H
#define cudasolver_H
/*----------------------------   cudasolver.h     ---------------------------*/

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <initializer_list>
#include <map>
#include <memory>
#include <stddef.h>
#include <string>

#include "check_cuda.h"
#include "cublas_v2.h"
#include "cudacsrmatrixinterface.h"
#include "cudavectorinterface.h"
#include "fmatrixblock.h"
#include "gascoigne.h"
#include "sparseblockmatrix.h"
#include "stdmultilevelsolver.h"
#include "stdsolver.h"
#include "vectorinterface.h"

namespace Gascoigne {
class CudaCSRMatrixInterface;
class MatrixInterface;

using CudaVectorAgent = GhostAgent<CudaVectorInterface>;

class CudaSolver : public StdSolver
{
protected:
  /// Cusparse handle
  cusparseHandle_t sparse_handle;
  cublasHandle_t blas_handle;

  /// Flag wheather or not the cuda is used. If false stdsolver is used.
  mutable bool oncuda = false;
  mutable CudaVectorAgent cva;

  IndexType nnodes = 0;

  mutable std::map<std::string, std::shared_ptr<CudaCSRMatrixInterface>>
    cuda_mat_agent;

public:
  CudaSolver();
  ~CudaSolver();

  void ActivateCuda(std::initializer_list<const Vector*>) const;
  void DeactivateCuda(std::initializer_list<Vector*> vectors) const;

  bool IsOnCuda() const { return oncuda; };

  void ReInitVector(Vector& dst) override;
  void DeleteVector(Vector& p) const override;

  CudaVectorInterface& InitCV(const Vector& u) const;
  CudaVectorInterface& GetCV(Vector& u) const;
  const CudaVectorInterface& GetCV(const Vector& u) const;
  virtual GlobalVector& GetGV(Vector& u) const;
  virtual const GlobalVector& GetGV(const Vector& u) const;

  std::shared_ptr<CudaCSRMatrixInterface> GetCudaMatrix(const Matrix& A) const;

  GlobalVector& CopyBack(Vector& u) const;

  void SetProblem(const ProblemDescriptorInterface& PDX) override;
  void AssembleMatrix(Matrix& A, Vector& u, double d) const override;

  void vmult(const Matrix& A,
             Vector& y,
             const Vector& x,
             double d) const override;
  void vmulteq(const Matrix& A,
               Vector& y,
               const Vector& x,
               double d) const override;

  void MatrixResidual(const Matrix& A,
                      Vector& gy,
                      const Vector& gx,
                      const Vector& gb) const override;
  void Jacobi(const Matrix& A, Vector& y) const override;
  void smooth(int niter,
              const Matrix& A,
              Vector& x,
              const Vector& y,
              Vector& h) const override;

  void Zero(Vector& dst) const override;

  void Equ(Vector& dst, double s, const Vector& src) const override;
  void Add(Vector& dst, double s, const Vector& src) const override;
  double Norm(const Vector& dst) const override;
  double ScalarProduct(const Vector& y, const Vector& x) const override;

  void HNZero(Vector& x) const override;
  void HNDistribute(Vector& x) const override;
  void HNAverage(Vector& x) const override;

  void residualgmres(const Matrix& A,
                     Vector& gy,
                     const Vector& gx,
                     const Vector& gb) const override;
  void SetBoundaryVectorZero(Vector& gf) const;
};
} // namespace Gascoigne

/*----------------------------   solver.h     ---------------------------*/
/* end of #ifndef cudasolver_H */
#endif
/*----------------------------   solver.h     ---------------------------*/
