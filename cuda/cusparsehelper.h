#ifndef CUSPARSEHELPER_H
#define CUSPARSEHELPER_H

#include <vector>

#include "columndiagstencil.h"
#include "fmatrixblock.h"

#include "check_cuda.h"

namespace Gascoigne {

/// Is not templated and runs the cuda part.
void
reorganise_host(size_t n,
                size_t n_entries,
                size_t n_row_max,
                size_t n_comp,
                MatrixEntryType* src_vals,
                size_t* src_rows,
                size_t* src_cols,
                MatrixEntryType* dest_vals,
                int32_t* dest_rows,
                int32_t* dest_cols);

template<int B>
void
reorganise(const std::vector<FMatrixBlock<B>>& data,
           const ColumnStencil* cs,
           MatrixEntryType* dest_vals,
           int32_t* dest_rows,
           int32_t* dest_cols)
{
  size_t n_comp = B;
  size_t n_entries = data.size();
  size_t n = cs->n();

  const IndexVector& start = cs->start();
  std::vector<int> start_int(start.begin(), start.end());

  const IndexVector& col = cs->col();
  std::vector<int> col_int(col.begin(), col.end());
  // size_t max_entries = 0;
  // for (size_t i = 0; i < n; ++i) {
  //   size_t div = start[i + 1] - start[i];
  //   if (div > max_entries) {
  //     max_entries = div;
  //   }
  // }

  MatrixEntryType* src_vals = nullptr;
  CHECK_CUDA(cudaMalloc(&src_vals, data.size() * sizeof(FMatrixBlock<B>)));
  CHECK_CUDA(cudaMemcpy(src_vals,
                        data.data(),
                        data.size() * sizeof(FMatrixBlock<B>),
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

  cusparseMatDescr_t descrA;
  CHECK_CUSPARSE(cusparseCreateMatDescr(&descrA));

  cusparseHandle_t sparse_handle;
  CHECK_CUSPARSE(cusparseCreate(&sparse_handle));
  CHECK_CUSPARSE(cusparseDbsr2csr(sparse_handle,
                                  CUSPARSE_DIRECTION_ROW,
                                  n,
                                  n,
                                  descrA,
                                  src_vals,
                                  src_rows,
                                  src_cols,
                                  n_comp,
                                  descrA,
                                  dest_vals,
                                  dest_rows,
                                  dest_cols));

  // reorganise_host(n,
  //                 n_entries,
  //                 n_comp,
  //                 max_entries,
  //                 src_vals,
  //                 src_rows,
  //                 src_cols,
  //                 dest_vals,
  //                 dest_rows,
  //                 dest_cols);

  CHECK_CUDA(cudaFree(src_vals));
  src_vals = nullptr;
  CHECK_CUDA(cudaFree(src_rows));
  src_rows = nullptr;
  CHECK_CUDA(cudaFree(src_cols));
  src_cols = nullptr;
}

// This is important for alignment issues.
#ifdef __MATRIX_SINGLE_PRECISION__
typedef uint32_t uint;
#else
typedef uint64_t uint;
#endif

/// Is not templated and runs the cuda part.
void
invert_host(size_t n_entries,
            size_t n_comp,
            MatrixEntryType* src_vals,
            IndexType* src_diag,
            MatrixEntryType* dest_vals,
            int32_t* dest_rows,
            int32_t* dest_cols);

template<int B>
void
invert(const std::vector<FMatrixBlock<B>>& data,
       const IndexVector& diag,
       MatrixEntryType* dest_vals,
       int32_t* dest_rows,
       int32_t* dest_cols)
{
  size_t n = diag.size();
  size_t n_comp = B;

  MatrixEntryType* src_vals;
  CHECK_CUDA(cudaMalloc(&src_vals, data.size() * sizeof(FMatrixBlock<B>)));
  CHECK_CUDA(cudaMemcpy(src_vals,
                        data.data(),
                        data.size() * sizeof(FMatrixBlock<B>),
                        cudaMemcpyHostToDevice));

  IndexType* src_diag;
  CHECK_CUDA(cudaMalloc(&src_diag, diag.size() * sizeof(IndexType)));
  CHECK_CUDA(cudaMemcpy(src_diag,
                        diag.data(),
                        diag.size() * sizeof(IndexType),
                        cudaMemcpyHostToDevice));

  invert_host(n, n_comp, src_vals, src_diag, dest_vals, dest_rows, dest_cols);

  CHECK_CUDA(cudaFree(src_vals));
  src_vals = nullptr;
  CHECK_CUDA(cudaFree(src_diag));
  src_diag = nullptr;
}

} // namespace Gascoigne
#endif