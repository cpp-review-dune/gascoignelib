/*----------------------------   cudainterface.h ---------------------------*/
/*      $Id:$                 */
#ifndef cudavectorinterface_H
#define cudavectorinterface_H
/*----------------------------   cudainterface.h ---------------------------*/

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <library_types.h>
#include <stddef.h>

#include "compvector.h"
#include "fmatrixblock.h"
#include "gascoigne.h"
#include "sparseblockmatrix.h"
#include "vectorinterface.h"

namespace Gascoigne {
template<class B>
class SparseBlockMatrix;
template<int N>
class FMatrixBlock;

/// Interface to access dense vectors on the device
struct CudaVectorInterface
{

  /// Device vector describtor for us in calculation
  cusparseDnVecDescr_t descr;

  /// Pointer to device memory of values
  double* valuesDevice = nullptr;

  /// Size of the vector
  size_t n;
  /// Size of single component
  size_t n_comp;

  /// Zero vector
  CudaVectorInterface();

  /// Zero vector
  CudaVectorInterface(size_t, size_t);

  CudaVectorInterface(const CudaVectorInterface&);
  CudaVectorInterface(CudaVectorInterface&&);

  CudaVectorInterface(const GlobalVector&);
  ~CudaVectorInterface();

  /// Copys data from rhs to this
  CudaVectorInterface& operator=(const GlobalVector&);
  /// Copys data from rhs to this
  CudaVectorInterface& operator=(const CudaVectorInterface&);

  void copy_back(GlobalVector&) const;

  /// Adds \par y to this.
  void add(cublasHandle_t, double, const CudaVectorInterface& y);

  /// Scals the vector by factor d
  void scal(cublasHandle_t, double d);

  /// Norm of the vector.
  double norm(cublasHandle_t) const;

  /// Zeros all data.
  void zero();
  /// Reserve a vector of lenght n * ncomp
  void reserve(size_t n, size_t n_comp);
  /// frees memory
  void free();

  operator bool() { return n * n_comp > 0; }
  friend std::ostream& operator<<(std::ostream&, const CudaVectorInterface&);
};

} // namespace Gascoigne

/*----------------------------   cudainterface.h ---------------------------*/
/* end of #ifndef cudavectorinterface_H */
#endif
/*----------------------------   cudainterface.h ---------------------------*/
