#include "cusparsehelper.h"

#include <vector>

#include <cuda_runtime.h>

#include "columnstencil.h"
#include "fmatrixblock.h"
#include "gascoigne.h"

#include "check_cuda.h"

namespace Gascoigne {

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