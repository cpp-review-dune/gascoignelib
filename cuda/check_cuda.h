#ifndef CHECK_CUDA_H
#define CHECK_CUDA_H

#include <cstdio> // cusparseSpMV
#include <cublas_v2.h>
#include <cuda_runtime.h> // cudaMalloc, cudaMemcpy, etc.
#include <cusparse.h>     // cusparseSpMV
#include <iostream>       // cusparseSpMV

#ifdef __unix__
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#endif

void inline print_stacktrace()
{
#ifdef __unix__
  void* array[10];
  int size;
  size = backtrace(array, 10);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
#endif
}

#define CHECK_CUDA(func)                                                       \
  {                                                                            \
    cudaError_t status = (func);                                               \
    if (status != cudaSuccess) {                                               \
      std::printf("CUDA API failed at %s:%d with error: %s (%d)\n",            \
                  __FILE__,                                                    \
                  __LINE__,                                                    \
                  cudaGetErrorString(status),                                  \
                  status);                                                     \
      print_stacktrace();                                                      \
      std::exit(1);                                                            \
    }                                                                          \
  }

#define CHECK_CUSPARSE(func)                                                   \
  {                                                                            \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
      std::printf("CUSPARSE API failed at %s:%d with error: %s (%d)\n",        \
                  __FILE__,                                                    \
                  __LINE__,                                                    \
                  cusparseGetErrorString(status),                              \
                  status);                                                     \
      print_stacktrace();                                                      \
      std::exit(1);                                                            \
    }                                                                          \
  }

#define CHECK_CUBLAS(func)                                                     \
  {                                                                            \
    cublasStatus_t status = (func);                                            \
    if (status != CUBLAS_STATUS_SUCCESS) {                                     \
      std::printf("CUBLAS API failed at %s:%d with error: (%d)\n",             \
                  __FILE__,                                                    \
                  __LINE__,                                                    \
                  status);                                                     \
      print_stacktrace();                                                      \
      std::exit(1);                                                            \
    }                                                                          \
  }

/**
 * Helpfull marker for debuging purpose
 */
#define CHECK                                                                  \
  {                                                                            \
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
  }

#endif // CHECK_CUDA_H
