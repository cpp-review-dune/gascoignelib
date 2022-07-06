#include "cusparsehelper.h"

#include <vector>

#include <cuda_runtime.h>

#include "columnstencil.h"
#include "fmatrixblock.h"
#include "gascoigne.h"

#include "check_cuda.h"

namespace Gascoigne {

__global__ void
reorganise_device(size_t n,
                  size_t n_entries,
                  size_t n_comp,
                  size_t n_row_max,
                  MatrixEntryType* src_vals,
                  size_t* src_rows,
                  size_t* src_cols,
                  MatrixEntryType* dest_vals,
                  int32_t* dest_rows,
                  int32_t* dest_cols)
{
  size_t id = threadIdx.x + blockIdx.x * blockDim.x;
  if (id >= n) {
    return;
  }

  // Copy Stuff to Shared Memory
  // 2*n_row_max*n_comp_sq * sizeof(MatrixEntryType) per thread
  extern __shared__ char sharedMemory[];
  size_t n_comp_sq = n_comp * n_comp;
  MatrixEntryType* shared_mat =
    (MatrixEntryType*)(sharedMemory) + 2 * threadIdx.x * n_row_max * n_comp_sq;
  uint* shared_cols = (uint*)shared_mat + n_row_max * n_comp_sq;

  size_t first_in_col = src_rows[id]; /// id of first entry in row
  size_t div =
    (src_rows[id + 1] - first_in_col) * n_comp; /// number of entries in row
  size_t dest_row = first_in_col * n_comp_sq;

  for (size_t i = 0; i < div * n_comp; ++i) {
    size_t matrix_id = i / (n_comp_sq);

    size_t in_matrix_id = i % (n_comp_sq);

    size_t rc = (in_matrix_id) / n_comp;
    size_t cc = (in_matrix_id) % n_comp;

    shared_mat[rc * div + matrix_id * n_comp + cc] = src_vals[dest_row + i];

    shared_cols[rc * div + matrix_id * n_comp + cc] =
      src_cols[first_in_col + i / n_comp_sq] * n_comp + cc;
  }

  for (size_t i = 0; i < div * n_comp; ++i) {
    dest_vals[dest_row + i] = shared_mat[i];
    dest_cols[dest_row + i] = shared_cols[i];
  }

  for (size_t i = 0; i < n_comp; ++i) {
    size_t org_row = (id * n_comp + i) / n_comp;
    size_t div = src_rows[org_row + 1] - src_rows[org_row];
    size_t org_div = (id * n_comp + i) % n_comp;
    dest_rows[(id * n_comp + i)] =
      src_rows[org_row] * n_comp_sq + org_div * div * n_comp;
  }

  if (id == 0) {
    dest_rows[n * n_comp] = src_rows[n] * n_comp_sq;
  }
}

void
reorganise_host(size_t n,
                size_t n_entries,
                size_t n_comp,
                size_t n_row_max,
                MatrixEntryType* src_vals,
                size_t* src_rows,
                size_t* src_cols,
                MatrixEntryType* dest_vals,
                int32_t* dest_rows,
                int32_t* dest_cols)
{
  if (n_entries == 0 || n == 0 || n_comp == 0) {
    return;
  }

  cudaDeviceProp props;
  CHECK_CUDA(cudaGetDeviceProperties(&props, 0));

  // size_t rows_per_block = 1024;

  // rows -> shared memory -> reorganise(row per thread) -> global memory

  size_t shared_per_row =
    n_comp * n_comp * n_row_max * sizeof(MatrixEntryType) * 2;
  size_t rows_per_block =
    min(static_cast<size_t>(floor((props.sharedMemPerBlock) /
                                  static_cast<double>(shared_per_row))),
        static_cast<size_t>(1024));
  size_t blocks = ceil((n_entries) / static_cast<double>(rows_per_block));

  reorganise_device<<<blocks,
                      rows_per_block,
                      shared_per_row * rows_per_block>>>(n,
                                                         n_entries,
                                                         n_comp,
                                                         n_row_max,
                                                         src_vals,
                                                         src_rows,
                                                         src_cols,
                                                         dest_vals,
                                                         dest_rows,
                                                         dest_cols);
  CHECK_CUDA(cudaDeviceSynchronize());
  CHECK_CUDA(cudaGetLastError());
}

__device__ void
copy_shared(MatrixEntryType* dest, MatrixEntryType* src, size_t size)
{
  for (size_t i = 0; i < size; ++i) {
    dest[i] = src[i];
  }
  __syncthreads();
}

__global__ void
invert_device(size_t n,
              size_t n_comp,
              size_t shared_per_mat,
              MatrixEntryType* src_vals,
              IndexType* src_diag,
              MatrixEntryType* dest_vals,
              int32_t* dest_rows,
              int32_t* dest_cols)
{
  extern __shared__ char sharedMemory[];
  size_t id = threadIdx.x + blockIdx.x * blockDim.x;
  if (id >= n) {
    return;
  }

  // assign shared memory
  MatrixEntryType* values =
    (MatrixEntryType*)(sharedMemory + threadIdx.x * shared_per_mat);
  MatrixEntryType* hv = (MatrixEntryType*)(values + n_comp * n_comp);
  uint* p = (uint*)(hv + n_comp);

  copy_shared(
    values, src_vals + n_comp * n_comp * src_diag[id], n_comp * n_comp);

  const size_t N = n_comp;

#define value(i, j) (values[i * N + j])

  int i, j, k, r;
  double max, hr;

  for (i = 0; i < N; i++)
    p[i] = i;

  for (j = 0; j < N; j++) {
    max = fabs(value(j, j));
    r = j;
    for (i = j + 1; i < N; i++) {
      if (fabs(value(i, j)) > max) {
        max = fabs(value(i, j));
        r = i;
      }
    }
    if (r > j) {
      for (k = 0; k < N; k++) {
        hr = value(j, k);
        value(j, k) = value(r, k);
        value(r, k) = hr;
      }
      i = p[j];
      p[j] = p[r];
      p[r] = i;
    }

    hr = 1. / value(j, j);
    value(j, j) = hr;
    for (k = 0; k < N; k++) {
      if (k == j)
        continue;
      for (i = 0; i < N; i++) {
        if (i == j)
          continue;
        value(i, k) -= value(i, j) * value(j, k) * hr;
      }
    }
    for (i = 0; i < N; i++) {
      value(i, j) *= hr;
      value(j, i) *= -hr;
    }
    value(j, j) = hr;
  }

  for (i = 0; i < N; i++) {
    for (k = 0; k < N; k++)
      hv[p[k]] = value(i, k);
    for (k = 0; k < N; k++)
      value(i, k) = hv[k];
  }

#undef value

  copy_shared(dest_vals + id * n_comp * n_comp, values, n_comp * n_comp);

  for (i = 0; i < n_comp; ++i) {
    for (j = 0; j < n_comp; ++j) {
      dest_cols[id * n_comp * n_comp + i * n_comp + j] = id * n_comp + j;
    }
    dest_rows[id * n_comp + i] = id * n_comp * n_comp + i * n_comp;
  }

  if (id == 0) {
    dest_rows[n * n_comp] = n * n_comp * n_comp;
  }
}

/// Is not templated and runs the cuda part.
void
invert_host(size_t n,
            size_t n_comp,
            MatrixEntryType* src_vals,
            IndexType* src_diag,
            MatrixEntryType* dest_vals,
            int32_t* dest_rows,
            int32_t* dest_cols)
{
  if (n == 0 || n_comp == 0) {
    return;
  }

  // shared memory of one thread needing the matrix and two helping vectors.
  size_t shared_per_mat = n_comp * n_comp * sizeof(MatrixEntryType) +
                          n_comp * sizeof(MatrixEntryType) +
                          n_comp * sizeof(uint);

  // Number of threads per block is limited by the memory or by 1024
  size_t threads_per_block =
    min(static_cast<size_t>(floor(32768 / static_cast<double>(shared_per_mat))),
        static_cast<size_t>(1024));
  size_t blocks =
    max(static_cast<size_t>(1),
        static_cast<size_t>(ceil(n / static_cast<double>(threads_per_block))));

  // Run kernal
  invert_device<<<blocks,
                  threads_per_block,
                  shared_per_mat * threads_per_block>>>(n,
                                                        n_comp,
                                                        shared_per_mat,
                                                        src_vals,
                                                        src_diag,
                                                        dest_vals,
                                                        dest_rows,
                                                        dest_cols);

  CHECK_CUDA(cudaDeviceSynchronize());
  CHECK_CUDA(cudaGetLastError());
}
} // namespace Gascoigne