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
class CudaCSRMatrixInterface
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

  CudaCSRMatrixInterface* dia_invers = nullptr;

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
  mutable bool bufferSet = false;
  mutable void* dBuffer = nullptr;

  /// Construct blank matrix, without any data

  CudaCSRMatrixInterface* get_inverse_diagonal(const MatrixEntryType* src_vals,
                                               const StencilInterface* stencil);

  template<int B>
  std::vector<FMatrixBlock<B>> patterToFMatrix(const std::vector<TimePattern>&);

  template<int B>
  void copy_data(const std::vector<FMatrixBlock<B>>& matrixdata,
                 const StencilInterface* stencil,
                 IndexType cols);

  template<int B>
  void copy_data(const SimpleMatrix&, const TimePattern*, IndexType);

  /**
   * Creates a inverted diagonal matrix of this one.
   */
  std::shared_ptr<CudaCSRMatrixInterface> get_inverse_diagonal(
    const CudaCSRMatrixInterface&);

  CudaCSRMatrixInterface(cusparseHandle_t sparse_handle,
                         IndexType n,
                         IndexType n_entries,
                         IndexType n_comp);

public:
  /// Construct matrix from MatrixInterface

  CudaCSRMatrixInterface(cusparseHandle_t, const MatrixInterface&);

  template<int B>
  CudaCSRMatrixInterface(cusparseHandle_t,
                         const std::vector<FMatrixBlock<B>>&,
                         const ColumnStencil*,
                         IndexType);

  CudaCSRMatrixInterface(cusparseHandle_t,
                         const std::vector<TimePattern>&,
                         const ColumnStencil*,
                         IndexType);

  CudaCSRMatrixInterface(cusparseHandle_t,
                         const SimpleMatrix&,
                         const TimePattern*,
                         IndexType);

  virtual ~CudaCSRMatrixInterface();

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

  void Jacobi(CudaVectorInterface& x) const;

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
  friend std::ostream& operator<<(std::ostream&, const CudaCSRMatrixInterface&);
};

// template<int B>
// std::ostream&
// operator<<(std::ostream& os, const CudaCSRMatrix<B>& dt);

} // namespace Gascoigne

// #include "cudacsrmatrixinterface.cc"

/*----------------------------   cudainterface.h ---------------------------*/
/* end of #ifndef cudainterface_H */
#endif
/*----------------------------   cudainterface.h ---------------------------*/
