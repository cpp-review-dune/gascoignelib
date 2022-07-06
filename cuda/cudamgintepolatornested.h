/*----------------------------   cudamgintepolatornested.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef cudamgintepolatornested_H
#define cudamgintepolatornested_H
/*----------------------------   cudamgintepolatornested.h
 * ---------------------------*/

#include <mginterpolatornested.h>
#include <scaledimatrixblock.h>
#include <simplematrix.h>

#include "cudacsrmatrixinterface.h"
#include "cudavectorinterface.h"
#include "gascoignemeshtransfer.h"

namespace Gascoigne {
class CudaMgInterpolatorNested : public MgInterpolatorNested
{
  /// Cusparse handle
  cusparseHandle_t sparse_handle;
  cublasHandle_t blas_handle;

  const GascoigneMeshTransfer* GT;

  mutable std::shared_ptr<SimpleMatrix> restrict;
  mutable std::shared_ptr<SimpleMatrix> prolongate;

  mutable size_t ncomp = 1;
  mutable std::string matrixtype;
  mutable std::shared_ptr<CudaCSRMatrixInterface> restrict_device = nullptr;
  mutable std::shared_ptr<CudaCSRMatrixInterface> prolongate_device = nullptr;

  void check_matrices(IndexType, IndexType H, IndexType h) const;

public:
  CudaMgInterpolatorNested();

  void BasicInit(const MeshTransferInterface* MT) override;
  void restrict_zero(GlobalVector&, const GlobalVector&) const override;
  void restrict_zero(CudaVectorInterface&, const CudaVectorInterface&) const;
  void prolongate_add(GlobalVector&, const GlobalVector&) const override;
  void prolongate_add(CudaVectorInterface&, const CudaVectorInterface&) const;
};

}
/*----------------------------   cudamgintepolatornested.h
 * ---------------------------*/
/* end of #ifndef cudamgintepolatornested */
#endif
/*----------------------------   cudamgintepolatornested.h
 * ---------------------------*/