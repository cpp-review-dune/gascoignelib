#include "cudacsrmatrixinterface.h"

#include <cassert>
#include <cstdint>
#include <cuda_runtime.h>
#include <ext/alloc_traits.h>
#include <memory>
#include <typeinfo>
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

CudaCSRMatrixInterface*
CudaCSRMatrixInterface::get_inverse_diagonal(const MatrixEntryType* src_vals,
                                             const StencilInterface* stencil)
{
  const ColumnDiagStencil* coldiag =
    dynamic_cast<const ColumnDiagStencil*>(stencil);
  if (!coldiag) {
    return nullptr;
  }
  const IndexVector& diag = coldiag->diag();
  IndexType n = diag.size();

  CudaCSRMatrixInterface* device_mat =
    new CudaCSRMatrixInterface(sparse_handle, n, n, _n_comp);

  // invert_direct(src_vals,
  //               diag,
  //               device_mat->csrValsDevice,
  //               device_mat->csrRowIdxsDevice,
  //               device_mat->csrColIdxsDevice,
  //               _n_comp);

  IndexType* src_diag;
  CHECK_CUDA(cudaMalloc(&src_diag, diag.size() * sizeof(IndexType)));
  CHECK_CUDA(cudaMemcpy(src_diag,
                        diag.data(),
                        diag.size() * sizeof(IndexType),
                        cudaMemcpyHostToDevice));

  invert_host(n,
              _n_comp,
              src_vals,
              src_diag,
              device_mat->csrValsDevice,
              device_mat->csrRowIdxsDevice,
              device_mat->csrColIdxsDevice);

  CHECK_CUDA(cudaFree(src_diag));
  src_diag = nullptr;

  return device_mat;
}

template<int B>
void
CudaCSRMatrixInterface::copy_data(
  const std::vector<FMatrixBlock<B>>& matrixdata,
  const StencilInterface* stencil,
  IndexType cols)
{
  GlobalTimer.count("---> cudamatrix mat");
  _n_comp = B;
  _n_entries = matrixdata.size();
  const ColumnStencil* CS = dynamic_cast<const ColumnStencil*>(stencil);
  if (!CS) {
    throw std::runtime_error("Wrong Stencil Type");
  }
  size_t n_rows = _n * _n_comp + 1;
  size_t n_values = _n_entries * _n_comp * _n_comp;
  CHECK_CUDA(cudaMalloc(&csrValsDevice, n_values * sizeof(MatrixEntryType)));
  CHECK_CUDA(cudaMalloc(&csrRowIdxsDevice, n_rows * sizeof(int32_t)));
  CHECK_CUDA(cudaMalloc(&csrColIdxsDevice, n_values * sizeof(int32_t)));

  // reorganise(matrixdata, CS, csrValsDevice, csrRowIdxsDevice,
  // csrColIdxsDevice);

  const IndexVector& start = CS->start();
  std::vector<int> start_int(start.begin(), start.end());

  const IndexVector& col = CS->col();
  std::vector<int> col_int(col.begin(), col.end());

  MatrixEntryType* src_vals = nullptr;
  CHECK_CUDA(
    cudaMalloc(&src_vals, matrixdata.size() * sizeof(FMatrixBlock<B>)));
  CHECK_CUDA(cudaMemcpy(src_vals,
                        matrixdata.data(),
                        matrixdata.size() * sizeof(FMatrixBlock<B>),
                        cudaMemcpyHostToDevice));

  int* src_rows = nullptr;
  CHECK_CUDA(cudaMalloc(&src_rows, start_int.size() * sizeof(int)));
  CHECK_CUDA(cudaMemcpy(src_rows,
                        start_int.data(),
                        start_int.size() * sizeof(int),
                        cudaMemcpyHostToDevice));

  int* src_cols = nullptr;
  CHECK_CUDA(cudaMalloc(&src_cols, col_int.size() * sizeof(int)));
  CHECK_CUDA(cudaMemcpy(src_cols,
                        col_int.data(),
                        col_int.size() * sizeof(int),
                        cudaMemcpyHostToDevice));

  dia_invers = get_inverse_diagonal(src_vals, stencil);

  cusparseMatDescr_t descrA;
  CHECK_CUSPARSE(cusparseCreateMatDescr(&descrA));

  cusparseHandle_t sparse_handle;
  CHECK_CUSPARSE(cusparseCreate(&sparse_handle));
  CHECK_CUSPARSE(cusparseDbsr2csr(sparse_handle,
                                  CUSPARSE_DIRECTION_ROW,
                                  static_cast<int>(CS->n()),
                                  static_cast<int>(CS->n()),
                                  descrA,
                                  src_vals,
                                  src_rows,
                                  src_cols,
                                  static_cast<int>(_n_comp),
                                  descrA,
                                  csrValsDevice,
                                  csrRowIdxsDevice,
                                  csrColIdxsDevice));

  CHECK_CUDA(cudaFree(src_vals));
  src_vals = nullptr;
  CHECK_CUDA(cudaFree(src_rows));
  src_rows = nullptr;
  CHECK_CUDA(cudaFree(src_cols));
  src_cols = nullptr;

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
}

template<int B>
std::vector<FMatrixBlock<B>>
CudaCSRMatrixInterface::patterToFMatrix(const std::vector<TimePattern>& pattern)
{
  std::vector<FMatrixBlock<B>> out(pattern.size());
  _n_comp = B;
  for (size_t i = 0; i < out.size(); ++i) {
    std::memcpy(out[i].data(),
                pattern[i].data(),
                _n_comp * _n_comp * sizeof(MatrixEntryType));
  }
  return out;
}

template<int B>
void
CudaCSRMatrixInterface::copy_data(const Gascoigne::SimpleMatrix& sbm,
                                  const TimePattern* pattern,
                                  IndexType cols)
{
  GlobalTimer.count("---> cudamatrix simplemt");
  const DoubleVector vals = sbm.GetValues();
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
    matrixdata, dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), cols);
}
CudaCSRMatrixInterface::CudaCSRMatrixInterface(cusparseHandle_t sparse_handle,
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
  cudaMemset(
    csrValsDevice, 0, static_cast<int32_t>(n_values * sizeof(MatrixEntryType)));
  CHECK_CUDA(cudaMalloc(&csrRowIdxsDevice, n_rows * sizeof(int32_t)));
  cudaMemset(
    csrRowIdxsDevice, 0, static_cast<int32_t>(n_rows * sizeof(int32_t)));
  CHECK_CUDA(cudaMalloc(&csrColIdxsDevice, n_values * sizeof(int32_t)));
  cudaMemset(
    csrColIdxsDevice, 0, static_cast<int32_t>(n_values * sizeof(int32_t)));

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

CudaCSRMatrixInterface::CudaCSRMatrixInterface(cusparseHandle_t sparse_handle,
                                               const MatrixInterface& mat)
  : sparse_handle(sparse_handle)
  , descr(nullptr)
  , _n(mat.GetStencil()->n())
{
  try {
    const auto& sbm =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<1>>&>(mat);
    copy_data(
      sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  try {
    const auto& sbm =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<2>>&>(mat);
    copy_data(
      sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  try {
    const auto& sbm =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<3>>&>(mat);
    copy_data(
      sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  try {
    const auto& sbm =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<4>>&>(mat);
    copy_data(
      sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  try {
    const auto& sbm =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<5>>&>(mat);
    copy_data(
      sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  try {
    const auto& sbm =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<6>>&>(mat);
    copy_data(
      sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  try {
    const auto& sbm =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<7>>&>(mat);
    copy_data(
      sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  try {
    const auto& sbm =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<8>>&>(mat);
    copy_data(
      sbm.mat(), dynamic_cast<const ColumnStencil*>(sbm.GetStencil()), _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  try {
    const auto& sm = dynamic_cast<const SimpleMatrix&>(mat);
    TimePattern tp(1, 1);
    copy_data<1>(sm, &tp, _n);
    // dia_invers = get_inverse_diagonal(mat.GetStencil());
    return;
  } catch (const std::bad_cast& e) {
    (void)e;
  }
  throw std::runtime_error(
    std::string("Can not construct CudaCSRMatrixInterface from ") +
    mat.GetName());
}

CudaCSRMatrixInterface::CudaCSRMatrixInterface(
  cusparseHandle_t sparse_handle,
  const Gascoigne::SimpleMatrix& sbm,
  const TimePattern* pattern,
  IndexType cols)
  : sparse_handle(sparse_handle)
  , descr(nullptr)
  , _n(sbm.GetStencil()->n())
{
  if (pattern->n() == 1) {
    copy_data<1>(sbm, pattern, cols);
    // dia_invers = get_inverse_diagonal(sbm->GetStencil());
  } else if (pattern->n() == 2) {
    copy_data<2>(sbm, pattern, cols);
    // dia_invers = get_inverse_diagonal(sbm->GetStencil());
  } else if (pattern->n() == 3) {
    copy_data<3>(sbm, pattern, cols);
    // dia_invers = get_inverse_diagonal(sbm->GetStencil());
  } else if (pattern->n() == 4) {
    copy_data<4>(sbm, pattern, cols);
    // dia_invers = get_inverse_diagonal(sbm->GetStencil());
  } else if (pattern->n() == 5) {
    copy_data<5>(sbm, pattern, cols);
    // dia_invers = get_inverse_diagonal(sbm->GetStencil());
  } else if (pattern->n() == 6) {
    copy_data<6>(sbm, pattern, cols);
    // dia_invers = get_inverse_diagonal(sbm->GetStencil());
  } else if (pattern->n() == 7) {
    copy_data<7>(sbm, pattern, cols);
    // dia_invers = get_inverse_diagonal(sbm->GetStencil());
  } else if (pattern->n() == 8) {
    copy_data<8>(sbm, pattern, cols);
    // dia_invers = get_inverse_diagonal(sbm->GetStencil());
  }
}

CudaCSRMatrixInterface::CudaCSRMatrixInterface(
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
    // dia_invers = get_inverse_diagonal(CS);
  } else if (pattern[0].n() == 2) {
    copy_data<2>(patterToFMatrix<2>(pattern), CS, cols);
    // dia_invers = get_inverse_diagonal(CS);
  } else if (pattern[0].n() == 3) {
    copy_data<3>(patterToFMatrix<3>(pattern), CS, cols);
    // dia_invers = get_inverse_diagonal(CS);
  } else if (pattern[0].n() == 4) {
    copy_data<4>(patterToFMatrix<4>(pattern), CS, cols);
    // dia_invers = get_inverse_diagonal(CS);
  } else if (pattern[0].n() == 5) {
    copy_data<5>(patterToFMatrix<5>(pattern), CS, cols);
    // dia_invers = get_inverse_diagonal(CS);
  } else if (pattern[0].n() == 6) {
    copy_data<6>(patterToFMatrix<6>(pattern), CS, cols);
    // dia_invers = get_inverse_diagonal(CS);
  } else if (pattern[0].n() == 7) {
    copy_data<7>(patterToFMatrix<7>(pattern), CS, cols);
    // dia_invers = get_inverse_diagonal(CS);
  } else if (pattern[0].n() == 8) {
    copy_data<8>(patterToFMatrix<8>(pattern), CS, cols);
    // dia_invers = get_inverse_diagonal(CS);
  }
}

template<int B>
CudaCSRMatrixInterface::CudaCSRMatrixInterface(
  cusparseHandle_t sparse_handle,
  const std::vector<FMatrixBlock<B>>& pattern,
  const ColumnStencil* CS,
  IndexType cols)
  : sparse_handle(sparse_handle)
  , descr(nullptr)
  , _n(cols)
{
  copy_data<B>(pattern, CS, cols);
  // dia_invers = get_inverse_diagonal(CS);
}

CudaCSRMatrixInterface::~CudaCSRMatrixInterface()
{

  GlobalTimer.count("---> cudamatrix dtor");
  CHECK_CUSPARSE(cusparseDestroySpMat(descr));

  CHECK_CUDA(cudaFree(csrRowIdxsDevice));
  csrRowIdxsDevice = nullptr;
  CHECK_CUDA(cudaFree(csrColIdxsDevice));
  csrColIdxsDevice = nullptr;
  CHECK_CUDA(cudaFree(csrValsDevice));
  csrValsDevice = nullptr;

  if (bufferSet) {
    CHECK_CUDA(cudaFree(dBuffer));
    dBuffer = nullptr;
  }

  if (dia_invers) {
    delete dia_invers;
    dia_invers = nullptr;
  }
}

/// y = d * A * x + b * y
void
CudaCSRMatrixInterface::vmult(CudaVectorInterface& cuda_gy,
                              const CudaVectorInterface& cuda_gx,
                              double d,
                              double b) const
{
  double alpha = d;
  double beta = b;

  cusparseOperation_t optype = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cusparseSpMVAlg_t algtype = CUSPARSE_MV_ALG_DEFAULT;
  if (!bufferSet) {
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
    bufferSet = true;
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

void
CudaCSRMatrixInterface::Jacobi(CudaVectorInterface& h) const
{
  if (!dia_invers) {
    throw std::runtime_error("No Diagonal Stencil");
  }
  CudaVectorInterface tmp(h.n, h.n_comp);
  dia_invers->vmult(tmp, h, 1, 0);
  h = tmp;
}

std::ostream&
operator<<(std::ostream& os, const CudaCSRMatrixInterface& dt)
{
  IndexType n = dt.n();
  IndexType n_comp = dt.n_comp();
  IndexType n_entry = dt.n_entries();

  os << "CudaCSRMatrixInterface Size (n n_comp n_entry n_value): "
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

template CudaCSRMatrixInterface::CudaCSRMatrixInterface<3>(
  cusparseHandle_t,
  const std::vector<FMatrixBlock<3>>&,
  const ColumnStencil*,
  IndexType);

template CudaCSRMatrixInterface::CudaCSRMatrixInterface<6>(
  cusparseHandle_t,
  const std::vector<FMatrixBlock<6>>&,
  const ColumnStencil*,
  IndexType);

} // namespace Gascoigne
