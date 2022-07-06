#include "cudabsrmatrixinterface.h"

#include <cassert>
#include <cstdint>
#include <cuda_runtime.h>
#include <ext/alloc_traits.h>
#include <memory>
#include <vector>

#include "check_cuda.h"
#include "columndiagstencil.h"
#include "columnstencil.h"
#include "cudavectorinterface.h"
#include "cusparse.h"
#include "driver_types.h"
#include "fmatrixblock.h"
#include "gascoigne.h"
#include "matrixinterface.h"
#include "stencilinterface.h"
#include "stopwatch.h"

namespace Gascoigne {
template<class B>
class SparseBlockMatrix;

extern Timer GlobalTimer;

template<int B>
void
CudaBSRMatrixInterface::copy_data(
  const std::vector<FMatrixBlock<B>>& matrixdata,
  const ColumnStencil* CS,
  IndexType cols)
{
  _n_comp = B;
  _n_entries = matrixdata.size();
  size_t n_rows = _n * _n_comp + 1;
  size_t n_values = _n_entries * _n_comp * _n_comp;
  CHECK_CUDA(cudaMalloc(&csrValsDevice, n_values * sizeof(MatrixEntryType)));
  CHECK_CUDA(cudaMalloc(&csrRowIdxsDevice, n_rows * sizeof(int32_t)));
  CHECK_CUDA(cudaMalloc(&csrColIdxsDevice, n_values * sizeof(int32_t)));

  if (cols == 0) {
    cols = _n;
  }

  CHECK_CUSPARSE(cusparseCreateCsr(&descr,
                                   static_cast<int64_t>(n_rows - 1),
                                   static_cast<int64_t>(cols * _n_comp),
                                   static_cast<int64_t>(n_values),
                                   csrRowIdxsDevice,
                                   csrColIdxsDevice,
                                   csrValsDevice,
                                   CUSPARSE_INDEX_32I,
                                   CUSPARSE_INDEX_32I,
                                   CUSPARSE_INDEX_BASE_ZERO,
                                   type));
};

template<int B>
std::vector<FMatrixBlock<B>>
CudaBSRMatrixInterface::patterToFMatrix(const std::vector<TimePattern>& pattern)
{
  std::vector<FMatrixBlock<B>> out(pattern.size());
  _n_comp = B;
  for (size_t i = 0; i < out.size(); ++i) {
    for (IndexType j = 0; j < _n_comp * _n_comp; ++j) {
      out[i][j] = pattern[i][j];
    }
  }
  return out;
}

template<int B>
void
CudaBSRMatrixInterface::copy_data(const SparseBlockMatrix<FMatrixBlock<B>>& sbm)
{
  _n_comp = B;
  GlobalTimer.count("---> cudamatrix mat");

  copy_data<B>(
    sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
};

template<int B>
void
CudaBSRMatrixInterface::copy_data(
  const std::shared_ptr<Gascoigne::SimpleMatrix>& sbm,
  const TimePattern* pattern,
  IndexType cols)
{
  GlobalTimer.count("---> cudamatrix simplemt");
  const DoubleVector vals = sbm->GetValues();
  _n_entries = vals.size();
  _n_comp = B;

  FMatrixBlock<B> identity;
  for (IndexType i = 0; i < _n_comp * _n_comp; ++i) {
    identity[i] = (*pattern)[i];
  }

  std::vector<FMatrixBlock<B>> matrixdata(_n_entries);

  for (IndexType i = 0; i < _n_entries; ++i) {
    matrixdata[i] = identity;
    matrixdata[i] *= vals[i];
  }

  copy_data<B>(
    matrixdata, dynamic_cast<const ColumnStencil*>(sbm->GetStencil()), cols);
};

/**
 * Creates inverted diagonal matrix for jacobian
 */

template<int B>
std::shared_ptr<CudaBSRMatrixInterface>
CudaBSRMatrixInterface::get_inverse_diagonal(
  cusparseHandle_t sparse_handle,
  const SparseBlockMatrix<FMatrixBlock<B>>& sbm)
{
  const ColumnDiagStencil* coldiag =
    dynamic_cast<const ColumnDiagStencil*>(sbm.GetStencil());
  const IndexVector& diag = coldiag->diag();
  IndexType n = diag.size();
  IndexType n_comp = B;

  auto device_mat =
    std::make_shared<CudaBSRMatrixInterface>(sparse_handle, n, n, n_comp);

  invert<B>(sbm.mat(),
            diag,
            device_mat->csrValsDevice,
            device_mat->csrRowIdxsDevice,
            device_mat->csrColIdxsDevice);

  return std::dynamic_pointer_cast<CudaBSRMatrixInterface>(device_mat);
};

CudaBSRMatrixInterface::CudaBSRMatrixInterface(cusparseHandle_t sparse_handle,
                                               IndexType n,
                                               IndexType n_entries,
                                               IndexType n_comp)
  : sparse_handle(sparse_handle)
  , descr(nullptr)
  , _n(n)
  , _n_entries(n_entries)
  , _n_comp(n_comp)
{
  GlobalTimer.count("---> cudamatrix size");

  size_t n_rows = n * _n_comp + 1;
  size_t n_values = n_entries * _n_comp * _n_comp;

  /* allocate GPU memory and copy the matrix and vectors into it */
  CHECK_CUDA(cudaMalloc(&csrValsDevice, n_values * sizeof(MatrixEntryType)));
  CHECK_CUDA(cudaMalloc(&csrRowIdxsDevice, n_rows * sizeof(int32_t)));
  CHECK_CUDA(cudaMalloc(&csrColIdxsDevice, n_values * sizeof(int32_t)));

  CHECK_CUSPARSE(cusparseCreateCsr(&descr,
                                   static_cast<int64_t>(n_rows - 1),
                                   static_cast<int64_t>(n_rows - 1),
                                   static_cast<int64_t>(n_values),
                                   csrRowIdxsDevice,
                                   csrColIdxsDevice,
                                   csrValsDevice,
                                   CUSPARSE_INDEX_32I,
                                   CUSPARSE_INDEX_32I,
                                   CUSPARSE_INDEX_BASE_ZERO,
                                   type));
}

CudaBSRMatrixInterface::CudaBSRMatrixInterface(cusparseHandle_t sparse_handle,
                                               const MatrixInterface& mat)
  : sparse_handle(sparse_handle)
  , descr(nullptr)
  , _n(mat.GetStencil()->n())
{
  try {
    copy_data<1>(dynamic_cast<const SparseBlockMatrix<FMatrixBlock<1>>&>(mat));
    return;
  } catch (std::bad_cast e) {
  }
  try {
    copy_data<2>(dynamic_cast<const SparseBlockMatrix<FMatrixBlock<2>>&>(mat));
    return;
  } catch (std::bad_cast e) {
  }
  try {
    copy_data<3>(dynamic_cast<const SparseBlockMatrix<FMatrixBlock<3>>&>(mat));
    return;
  } catch (std::bad_cast e) {
  }
  try {
    copy_data<4>(dynamic_cast<const SparseBlockMatrix<FMatrixBlock<4>>&>(mat));
    return;
  } catch (std::bad_cast e) {
  }
  try {
    copy_data<5>(dynamic_cast<const SparseBlockMatrix<FMatrixBlock<5>>&>(mat));
    return;
  } catch (std::bad_cast e) {
  }
  try {
    copy_data<6>(dynamic_cast<const SparseBlockMatrix<FMatrixBlock<6>>&>(mat));
    return;
  } catch (std::bad_cast e) {
  }
  try {
    copy_data<7>(dynamic_cast<const SparseBlockMatrix<FMatrixBlock<7>>&>(mat));
    return;
  } catch (std::bad_cast e) {
  }
  try {
    copy_data<8>(dynamic_cast<const SparseBlockMatrix<FMatrixBlock<8>>&>(mat));
    return;
  } catch (std::bad_cast e) {
  }
}

CudaBSRMatrixInterface::CudaBSRMatrixInterface(
  cusparseHandle_t sparse_handle,
  const std::shared_ptr<Gascoigne::SimpleMatrix>& sbm,
  const TimePattern* pattern,
  IndexType cols)
  : sparse_handle(sparse_handle)
  , descr(nullptr)
  , _n(sbm->GetStencil()->n())
{
  if (pattern->n() == 1) {
    copy_data<1>(sbm, pattern, cols);
  } else if (pattern->n() == 2) {
    copy_data<2>(sbm, pattern, cols);
  } else if (pattern->n() == 3) {
    copy_data<3>(sbm, pattern, cols);
  } else if (pattern->n() == 4) {
    copy_data<4>(sbm, pattern, cols);
  } else if (pattern->n() == 5) {
    copy_data<5>(sbm, pattern, cols);
  } else if (pattern->n() == 6) {
    copy_data<6>(sbm, pattern, cols);
  } else if (pattern->n() == 7) {
    copy_data<7>(sbm, pattern, cols);
  } else if (pattern->n() == 8) {
    copy_data<8>(sbm, pattern, cols);
  }
}

CudaBSRMatrixInterface::CudaBSRMatrixInterface(
  cusparseHandle_t sparse_handle,
  const std::vector<TimePattern>& pattern,
  const ColumnStencil* CS,
  IndexType cols)
  : sparse_handle(sparse_handle)
  , descr(nullptr)
  , _n(cols)
{
  if (pattern[0].n() == 1) {
    copy_data<1>(patterToFMatrix<1>(pattern), CS, cols);
  } else if (pattern[0].n() == 2) {
    copy_data<2>(patterToFMatrix<2>(pattern), CS, cols);
  } else if (pattern[0].n() == 3) {
    copy_data<3>(patterToFMatrix<3>(pattern), CS, cols);
  } else if (pattern[0].n() == 4) {
    copy_data<4>(patterToFMatrix<4>(pattern), CS, cols);
  } else if (pattern[0].n() == 5) {
    copy_data<5>(patterToFMatrix<5>(pattern), CS, cols);
  } else if (pattern[0].n() == 6) {
    copy_data<6>(patterToFMatrix<6>(pattern), CS, cols);
  } else if (pattern[0].n() == 7) {
    copy_data<7>(patterToFMatrix<7>(pattern), CS, cols);
  } else if (pattern[0].n() == 8) {
    copy_data<8>(patterToFMatrix<8>(pattern), CS, cols);
  }
}

/**
 * Creates inverted diagonal matrix for jacobian
 */
std::shared_ptr<CudaBSRMatrixInterface>
CudaBSRMatrixInterface::get_inverse_diagonal(cusparseHandle_t sparse_handle,
                                             const MatrixInterface& mat)
{
  try {
    return get_inverse_diagonal<1>(
      sparse_handle,
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<1>>&>(mat));
  } catch (std::bad_cast e) {
  }
  try {
    return get_inverse_diagonal<2>(
      sparse_handle,
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<2>>&>(mat));
  } catch (std::bad_cast e) {
  }
  try {
    return get_inverse_diagonal<3>(
      sparse_handle,
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<3>>&>(mat));
  } catch (std::bad_cast e) {
  }
  try {
    return get_inverse_diagonal<4>(
      sparse_handle,
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<4>>&>(mat));
  } catch (std::bad_cast e) {
  }
  try {
    return get_inverse_diagonal<5>(
      sparse_handle,
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<5>>&>(mat));
  } catch (std::bad_cast e) {
  }
  try {
    return get_inverse_diagonal<6>(
      sparse_handle,
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<6>>&>(mat));
  } catch (std::bad_cast e) {
  }
  try {
    return get_inverse_diagonal<7>(
      sparse_handle,
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<7>>&>(mat));
  } catch (std::bad_cast e) {
  }
  try {
    return get_inverse_diagonal<8>(
      sparse_handle,
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<8>>&>(mat));
  } catch (std::bad_cast e) {
  }
  abort();
}

CudaBSRMatrixInterface::~CudaBSRMatrixInterface()
{
  GlobalTimer.count("---> cudamatrix dtor");
  CHECK_CUSPARSE(cusparseDestroySpMat(descr));

  CHECK_CUDA(cudaFree(csrRowIdxsDevice));
  csrRowIdxsDevice = nullptr;
  CHECK_CUDA(cudaFree(csrColIdxsDevice));
  csrColIdxsDevice = nullptr;
  CHECK_CUDA(cudaFree(csrValsDevice));
  csrValsDevice = nullptr;

  if (dBuffer) {
    CHECK_CUDA(cudaFree(dBuffer));
    dBuffer = nullptr;
  }
}

/// y = d * A * x + b * y
void
CudaBSRMatrixInterface::vmult(CudaVectorInterface& cuda_gy,
                              const CudaVectorInterface& cuda_gx,
                              double d,
                              double b) const
{
  double alpha = d;
  double beta = b;

  cusparseOperation_t optype = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cusparseSpMVAlg_t algtype = CUSPARSE_MV_ALG_DEFAULT;
  if (!dBuffer) {
    size_t bufferSize = 0;
    CHECK_CUSPARSE(cusparseSpMV_bufferSize(sparse_handle,
                                           optype,
                                           &alpha,
                                           descr,
                                           cuda_gx.descr,
                                           &beta,
                                           cuda_gy.descr,
                                           type,
                                           algtype,
                                           &bufferSize));

    CHECK_CUDA(cudaMalloc(&dBuffer, bufferSize));
  }

  // cuda_gy = alpha * descr * cuda_gx + beta * cuda_gy
  CHECK_CUSPARSE(cusparseSpMV(sparse_handle,
                              optype,
                              &alpha,
                              descr,
                              cuda_gx.descr,
                              &beta,
                              cuda_gy.descr,
                              type,
                              algtype,
                              dBuffer));
  CHECK_CUDA(cudaDeviceSynchronize());
}

// template<int B>
// void
// CudaBSRMatrixInterface::operator=(const std::vector<FMatrixBlock<B>>&
// matrixdata)
// {
//   size_t n_values = _n_entries * _n_comp * _n_comp;

//   std::vector<MatrixEntryType> csrVals(n_values);

//   for (size_t i = 0; i < _n; ++i) { // schleife ueber alle Zeilen
//     FMatrixBlock<B> b = matrixdata[i];
//     for (size_t rc = 0; rc < _n_comp; ++rc) {   // component-row //
//       for (size_t cc = 0; cc < _n_comp; ++cc) { // component-row //
//         csrVals[i * _n_comp * _n_comp + rc * _n_comp + cc] =
//           b(static_cast<int>(rc), static_cast<int>(cc));
//       }
//     }
//   }
//   CHECK_CUDA(cudaMemcpy(csrValsDevice,
//                         csrVals.data(),
//                         n_values * sizeof(csrVals[0]),
//                         cudaMemcpyHostToDevice));
// }

std::ostream&
operator<<(std::ostream& os, const CudaBSRMatrixInterface& dt)
{
  IndexType n = dt.n();
  IndexType n_comp = dt.n_comp();
  IndexType n_entry = dt.n_entries();

  os << "CudaBSRMatrixInterface Size (n n_comp n_entry n_value): "
     << std::to_string(n) << " " << std::to_string(n_comp) << " "
     << std::to_string(n_entry) << " "
     << std::to_string(n_entry * n_comp * n_comp) << std::endl;

  if (n == 0 || n_comp == 0 || n_entry == 0) {
    os << "Empty matrix!" << std::endl;
    return os;
  }

  std::vector<MatrixEntryType> diag_check(n_entry * n_comp * n_comp);
  CHECK_CUDA(
    cudaMemcpy(diag_check.data(),
               dt.csrValsDevice,
               static_cast<int32_t>(diag_check.size() * sizeof(double)),
               cudaMemcpyDeviceToHost))
  std::vector<int32_t> diag_row_check(n * n_comp + 1);
  CHECK_CUDA(
    cudaMemcpy(diag_row_check.data(),
               dt.csrRowIdxsDevice,
               static_cast<int32_t>(diag_row_check.size() * sizeof(int32_t)),
               cudaMemcpyDeviceToHost))
  std::vector<int32_t> diag_cols_check(n_entry * n_comp * n_comp);
  CHECK_CUDA(
    cudaMemcpy(diag_cols_check.data(),
               dt.csrColIdxsDevice,
               static_cast<int32_t>(diag_cols_check.size() * sizeof(int32_t)),
               cudaMemcpyDeviceToHost))

  os << "Row Starts: ";
  auto row_iter = diag_row_check.begin();
  for (size_t i = 0; i < 20; ++i) {
    os << std::setw(9) << std::right << (row_iter[i]) << " ";
  }
  os << std::endl;
  for (size_t i = 0; i < 20; ++i) {
    os << std::setw(9) << std::right << diag_cols_check[i] << " "
       << std::setw(9) << std::right << diag_check[i] << std::endl;
  }
  return os;
}
} // namespace Gascoigne
