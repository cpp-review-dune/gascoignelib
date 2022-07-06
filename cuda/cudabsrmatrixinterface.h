/*----------------------------   cudainterface.h ---------------------------*/
/*      $Id:$                 */
#ifndef cudainterface_H
#define cudainterface_H
/*----------------------------   cudainterface.h ---------------------------*/

#include <cuda_runtime.h>
#include <cusparse.h>
#include <iostream>
#include <library_types.h>
#include <memory>
#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "compvector.h"
#include "cudavectorinterface.h"
#include "cusparsehelper.h"
#include "fmatrixblock.h"
#include "gascoigne.h"
#include "simplematrix.h"
#include "sparseblockmatrix.h"
#include "stopwatch.h"

namespace Gascoigne {
class MatrixInterface;
template<int N>
class FMatrixBlock;
extern Timer GlobalTimer;

/** Interface to access csr matrix on the device
 *
 * n_comp     := B
 * n          := the size of the matrix in FMatrixBlock's
 * n_entries  := the number of FMatrixBlock's
 *
 * n_values   := n_entries * n_comp * n_comp
 * n_rows     := n * n_comp + 1
 */
class CudaBSRMatrixInterface
{
private:
  /// cusparse handle needed for calculations.
  cusparseHandle_t sparse_handle;

  /// Cusparse descriptor for matrix
  cusparseSpMatDescr_t descr;

  /// Cuda type, may be CUDA_R_64F or CUDA_R_32F
#ifdef __MATRIX_SINGLE_PRECISION__
  const cudaDataType type = CUDA_R_32F;
#else
  const cudaDataType type = CUDA_R_64F;
#endif

  size_t _n;
  size_t _n_entries;
  size_t _n_comp;

  /** Pointer to Device Memory for values.
   * csrColIdxsDevice.size() == n_entries * n_comp * n_comp
   */
  MatrixEntryType* csrValsDevice = nullptr;

  /** Pointer to Device Memory for row indexes.
   * csrRowIdxsDevice.size() == n_rows == n * n_comp +1
   */
  int32_t* csrRowIdxsDevice = nullptr;

  /** Pointer to Device Memory for column indexes.
   * csrColIdxsDevice.size() == n_entries * n_comp * n_comp
   */
  int32_t* csrColIdxsDevice = nullptr;

  /**
   * Buffer for vmult on device
   */
  mutable void* dBuffer = nullptr;

  /// Construct blank matrix, without any data

  template<int B>
  static std::shared_ptr<CudaBSRMatrixInterface> get_inverse_diagonal(
    cusparseHandle_t,
    const SparseBlockMatrix<FMatrixBlock<B>>&);

  template<int B>
  std::vector<FMatrixBlock<B>> patterToFMatrix(const std::vector<TimePattern>&);

  template<int B>
  void copy_data(const std::vector<FMatrixBlock<B>>&,
                 const ColumnStencil*,
                 IndexType);

  template<int B>
  void copy_data(const SparseBlockMatrix<FMatrixBlock<B>>&);

  template<int B>
  void copy_data(const std::shared_ptr<Gascoigne::SimpleMatrix>&,
                 const TimePattern*,
                 IndexType);

public:
  CudaBSRMatrixInterface(cusparseHandle_t sparse_handle,
                         IndexType n,
                         IndexType n_entries,
                         IndexType n_comp);
  /// Construct matrix from MatrixInterface

  CudaBSRMatrixInterface(cusparseHandle_t, const MatrixInterface&);
  CudaBSRMatrixInterface(cusparseHandle_t,
                         const std::vector<TimePattern>&,
                         const ColumnStencil*,
                         IndexType);
  CudaBSRMatrixInterface(cusparseHandle_t,
                         const std::shared_ptr<Gascoigne::SimpleMatrix>&,
                         const TimePattern*,
                         IndexType);

  // CudaBSRMatrixInterface(cusparseHandle_t sparse_handle,
  //               const std::vector<MatrixEntryType> &csrVals,
  //               const std::vector<int32_t> &csrRowIdxs,
  //               const std::vector<int32_t> &csrColIdxs);
  virtual ~CudaBSRMatrixInterface();

  /**
   * Creates a inverted diagonal matrix of this one.
   */
  static std::shared_ptr<CudaBSRMatrixInterface> get_inverse_diagonal(
    cusparseHandle_t,
    const MatrixInterface&);

  size_t n() const { return _n; }
  size_t n_comp() const { return _n_comp; };
  size_t n_entries() const { return _n_entries; };

  /**
   * vmult function from superclass with calculations on cuda device
   *  y = d * A * x + b * y
   */
  void vmult(CudaVectorInterface& y,
             const CudaVectorInterface& x,
             double d,
             double b) const;

  /**
   * Operator to overwrite the values;
   */
  // void operator=(const std::vector<FMatrixBlock<B>>& matrixdata);

  // private:
  //   void copy_to_device(const std::vector<MatrixEntryType> &csrVals,
  //                       const std::vector<int32_t> &csrRowIdxs,
  //                       const std::vector<int32_t> &csrColIdxs);
  //   void copy_to_device(size_t n_rows, size_t n_values,
  //                       const MatrixEntryType *csrVals, const int32_t
  //                       *csrRowIdxs, const int32_t *csrColIdxs);
  friend std::ostream& operator<<(std::ostream&, const CudaBSRMatrixInterface&);
};

// template<int B>
// std::ostream&
// operator<<(std::ostream& os, const CudaCSRMatrix<B>& dt);

} // namespace Gascoigne

// #include "CudaBSRMatrixInterface.cc"

/*----------------------------   cudainterface.h ---------------------------*/
/* end of #ifndef cudainterface_H */
#endif
/*----------------------------   cudainterface.h ---------------------------*/
