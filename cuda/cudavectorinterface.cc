#include "cudavectorinterface.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <library_types.h>
#include <stdint.h>

#include "check_cuda.h"
#include "driver_types.h"
#include "stopwatch.h"

namespace Gascoigne {
extern Timer GlobalTimer;

CudaVectorInterface::CudaVectorInterface()
  : CudaVectorInterface(0, 1)
{
  GlobalTimer.count("---> cuvec empty");
}

CudaVectorInterface::CudaVectorInterface(size_t n, size_t n_comp)
  : descr()
  , valuesDevice(nullptr)
  , n(0)
  , n_comp(1)
{
  GlobalTimer.count("---> cuvec size");
  if (n * n_comp > 0) {
    reserve(n, n_comp);
    zero();
  }
}

CudaVectorInterface::CudaVectorInterface(const CudaVectorInterface& src)
  : descr()
  , valuesDevice(nullptr)
  , n(0)
  , n_comp(1)
{
  GlobalTimer.count("---> cuvec cuvec");
  reserve(src.n, src.n_comp);

  // Copy memory to device
  CHECK_CUDA(cudaMemcpy(valuesDevice,
                        src.valuesDevice,
                        static_cast<int32_t>(n * n_comp * sizeof(double)),
                        cudaMemcpyDeviceToDevice));
}

CudaVectorInterface::CudaVectorInterface(CudaVectorInterface&& src)
  : descr(src.descr)
  , valuesDevice(src.valuesDevice)
  , n(src.n)
  , n_comp(src.n_comp)
{
  GlobalTimer.count("---> cuvec cuvec&&");
  // src.valuesDevice = nullptr;
}

CudaVectorInterface::CudaVectorInterface(const GlobalVector& host_vector)
  : descr()
  , valuesDevice(nullptr)
  , n(0)
  , n_comp(0)
{
  GlobalTimer.count("---> cuvec glob");
  reserve(host_vector.n(), host_vector.ncomp());
  // Copy memory to device
  CHECK_CUDA(cudaMemcpy(valuesDevice,
                        host_vector.data(),
                        static_cast<int32_t>(n * n_comp * sizeof(double)),
                        cudaMemcpyHostToDevice));
}

CudaVectorInterface&
CudaVectorInterface::operator=(const GlobalVector& rhs)
{
  GlobalTimer.count("---> cuvec =glob");

  reserve(rhs.n(), rhs.ncomp());

  CHECK_CUDA(cudaMemcpy(valuesDevice,
                        rhs.data(),
                        static_cast<int32_t>(n * n_comp * sizeof(double)),
                        cudaMemcpyHostToDevice));
  return *this;
}

/// Copys data from rhs to this
CudaVectorInterface&
CudaVectorInterface::operator=(const CudaVectorInterface& rhs)
{
  GlobalTimer.count("---> cuvec =cuvec");
  reserve(rhs.n, rhs.n_comp);

  // Copy memory to device
  CHECK_CUDA(cudaMemcpy(valuesDevice,
                        rhs.valuesDevice,
                        static_cast<int32_t>(n * n_comp * sizeof(double)),
                        cudaMemcpyDeviceToDevice));
  return *this;
}

CudaVectorInterface::~CudaVectorInterface()
{
  free();
}

void
CudaVectorInterface::copy_back(GlobalVector& host_vector) const
{
  CHECK_CUDA(cudaMemcpy(host_vector.data(),
                        valuesDevice,
                        static_cast<int32_t>(n * n_comp * sizeof(double)),
                        cudaMemcpyDeviceToHost))
}

void
CudaVectorInterface::zero()
{
  if (valuesDevice != nullptr) {
    cudaMemset(
      valuesDevice, 0, static_cast<int32_t>(n * n_comp * sizeof(double)));
  }
}

void
CudaVectorInterface::free()
{
  // if (valuesDevice != nullptr) {
  CHECK_CUSPARSE(cusparseDestroyDnVec(descr));
  CHECK_CUDA(cudaFree(valuesDevice));
  GlobalTimer.count("---> cuvec frees");
  valuesDevice = nullptr;
  // } else {
  //   GlobalTimer.count("---> cuvec not frees");
  // }
}

void
CudaVectorInterface::reserve(size_t n, size_t n_comp)
{
  if (n == 0 || n_comp == 0) {
    free();
  } else if ((this->n * this->n_comp) != (n * n_comp)) {
    free();

    CHECK_CUDA(cudaMalloc(&valuesDevice, n * n_comp * sizeof(double)));
    GlobalTimer.count("---> cuvec malloc");

    CHECK_CUSPARSE(cusparseCreateDnVec(
      &descr, n * n_comp, valuesDevice, cudaDataType::CUDA_R_64F));
  }
  this->n = n;
  this->n_comp = n_comp;
}

/**
 * this = this + d * v
 */
void
CudaVectorInterface::add(cublasHandle_t handle,
                         double d,
                         const CudaVectorInterface& v)
{
  CHECK_CUBLAS(cublasDaxpy(handle,
                           static_cast<int>(n * n_comp),
                           &d,
                           v.valuesDevice,
                           1,
                           valuesDevice,
                           1));
}

double
CudaVectorInterface::norm(cublasHandle_t handle) const
{
  double result = 0;
  CHECK_CUBLAS(cublasDnrm2(
    handle, static_cast<int>(n * n_comp), valuesDevice, 1, &result));
  return result;
}

void
CudaVectorInterface::scal(cublasHandle_t handle, double d)
{
  CHECK_CUBLAS(
    cublasDscal(handle, static_cast<int>(n * n_comp), &d, valuesDevice, 1));
}

std::ostream&
operator<<(std::ostream& os, const CudaVectorInterface& vi)
{
  IndexType n = vi.n;
  IndexType n_comp = vi.n_comp;

  os << "CudaVectorInterface Size (n n_comp): " << std::to_string(n) << " "
     << std::to_string(n_comp) << std::endl;

  if (n == 0 || n_comp == 0) {
    os << "Empty Vector!";
    return os;
  }

  std::vector<double> value_check(n * n_comp);
  CHECK_CUDA(
    cudaMemcpy(value_check.data(),
               vi.valuesDevice,
               static_cast<int32_t>(value_check.size() * sizeof(double)),
               cudaMemcpyDeviceToHost));

  os << "Values: ";
  for (size_t i = 0; i < 20; ++i) {
    os << std::setw(9) << std::right << (value_check[i]) << " ";
  }
  return os;
}

} // namespace Gascoigne
