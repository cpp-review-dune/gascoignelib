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

  IndexType ncomp = 1;
  IndexType nnodes = 0;
  std::string matrixtype;

  mutable std::shared_ptr<CudaCSRMatrixInterface> mat = nullptr;
  mutable std::shared_ptr<CudaCSRMatrixInterface> dia_invers = nullptr;

  std::shared_ptr<SimpleMatrix> hn_zero;
  std::shared_ptr<SimpleMatrix> hn_average;
  std::shared_ptr<SimpleMatrix> hn_distribute;
  std::shared_ptr<CudaCSRMatrixInterface> hn_zero_device = nullptr;
  std::shared_ptr<CudaCSRMatrixInterface> hn_average_device = nullptr;
  std::shared_ptr<CudaCSRMatrixInterface> hn_distribute_device = nullptr;

  std::shared_ptr<CudaCSRMatrixInterface> dirichlet_zeros;

public:
  CudaSolver();
  ~CudaSolver();

  void ActivateCuda(std::initializer_list<const Vector*>) const;
  void DeactivateCuda(std::initializer_list<Vector*> vectors) const;

  bool IsOnCuda() const { return oncuda; };

  CudaVectorInterface& InitCV(const Vector& u) const;
  CudaVectorInterface& GetCV(const Vector& u) const;
  GlobalVector& CopyBack(Vector& u) const;
  void DeleteVector(Vector& p) const override;

  void BasicInit(const ParamFile& paramfile, const int dimension) override;
  MatrixInterface* NewMatrix(int ncomp, const std::string& matrixtype) override;
  void NewMesh(const GascoigneMesh* mp) override;
  void SetProblem(const ProblemDescriptorInterface& PDX) override;
  void AssembleMatrix(Matrix& A, const Vector& u, double d) const override;

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
  void smooth(int niter,
              const Matrix& A,
              Vector& x,
              const Vector& y,
              Vector& h) const override;

  void Zero(Vector& dst) const override;

  void Equ(Vector& dst, double s, const Vector& src) const override;
  void Add(Vector& dst, double s, const Vector& src) const override;
  double Norm(const Vector& dst) const override;

  void HNZero(const Vector& x) const override;
  void HNDistribute(Vector& x) const override;
  void HNAverage(const Vector& x) const override;
  void SetBoundaryVectorZero(Vector& gf) const;
};
} // namespace Gascoigne

/*----------------------------   solver.h     ---------------------------*/
/* end of #ifndef cudasolver_H */
#endif
/*----------------------------   solver.h     ---------------------------*/
