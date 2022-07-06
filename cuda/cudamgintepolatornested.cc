
#include "cudamgintepolatornested.h"

#include <algorithm>
#include <array>
#include <map>

#include "check_cuda.h"
#include "stopwatch.h"

namespace Gascoigne {
extern Timer GlobalTimer;

CudaMgInterpolatorNested::CudaMgInterpolatorNested()
  : MgInterpolatorNested()
{
  CHECK_CUSPARSE(cusparseCreate(&sparse_handle));
  CHECK_CUBLAS(cublasCreate(&blas_handle));
}

void
CudaMgInterpolatorNested::BasicInit(const MeshTransferInterface* MT)
{
  GT = dynamic_cast<const GascoigneMeshTransfer*>(MT);
}

void
CudaMgInterpolatorNested::check_matrices(IndexType n_comp,
                                         IndexType H,
                                         IndexType h) const
{
  if (prolongate == nullptr) {
    IntVector c2f = GT->GetC2f();
    std::map<int, std::array<int, 2>> zweier = GT->GetZweier();
    std::map<int, std::array<int, 4>> vierer = GT->GetVierer();
    std::map<int, std::array<int, 8>> achter = GT->GetAchter();

    SparseStructure saRestrict;
    saRestrict.build_begin(H);
    SparseStructure saProlongate;
    saProlongate.build_begin(h);

    for (int i = 0; i < c2f.size(); i++) {
      saRestrict.build_add(i, c2f[i]);
      saProlongate.build_add(c2f[i], i);
    }
    for (auto p = zweier.begin(); p != zweier.end(); p++) {
      int il = p->first;
      for (const auto& n : p->second) {
        saRestrict.build_add(n, il);
        saProlongate.build_add(il, n);
      }
    }
    for (auto p = vierer.begin(); p != vierer.end(); p++) {
      int il = p->first;
      for (const auto& n : p->second) {
        saRestrict.build_add(n, il);
        saProlongate.build_add(il, n);
      }
    }
    for (auto p = achter.begin(); p != achter.end(); p++) {
      int il = p->first;
      for (const auto& n : p->second) {
        saRestrict.build_add(n, il);
        saProlongate.build_add(il, n);
      }
    }
    saRestrict.build_end();
    saProlongate.build_end();

    restrict = std::make_shared<SimpleMatrix>();
    prolongate = std::make_shared<SimpleMatrix>();

    restrict->ReInit(&saRestrict);
    prolongate->ReInit(&saProlongate);

    for (int i = 0; i < c2f.size(); i++) {
      restrict->GetValue(i, c2f[i]) = 1;
      prolongate->GetValue(c2f[i], i) = 1;
    }
    for (auto p = zweier.begin(); p != zweier.end(); p++) {
      int il = p->first;
      for (const auto& n : p->second) {
        restrict->GetValue(n, il) = 0.5;
        prolongate->GetValue(il, n) = 0.5;
      }
    }
    for (auto p = vierer.begin(); p != vierer.end(); p++) {
      int il = p->first;
      for (const auto& n : p->second) {
        restrict->GetValue(n, il) = 0.25;
        prolongate->GetValue(il, n) = 0.25;
      }
    }
    for (auto p = achter.begin(); p != achter.end(); p++) {
      int il = p->first;
      for (const auto& n : p->second) {
        restrict->GetValue(n, il) = 0.125;
        prolongate->GetValue(il, n) = 0.125;
      }
    }

    TimePattern pattern(n_comp);
    pattern.identity();
    restrict_device = std::make_shared<CudaCSRMatrixInterface>(
      sparse_handle, restrict, &pattern, prolongate->GetStencil()->n());
    prolongate_device = std::make_shared<CudaCSRMatrixInterface>(
      sparse_handle, prolongate, &pattern, restrict->GetStencil()->n());
  }
}

// x_h -> x_H
// x_h > x_H
// x_H = R * x_h
// R \in \R^{H x h}
void
CudaMgInterpolatorNested::restrict_zero(GlobalVector& x_H,
                                        const GlobalVector& x_h) const
{
  GlobalTimer.count("---> restrict host");
  TimePattern pattern(x_H.ncomp());
  pattern.identity();
  check_matrices(x_H.ncomp(), x_H.n(), x_h.n());
  x_H.zero();
  restrict->vmult_time(x_H, x_h, pattern, 1);
}

// x_h -> x_H
// x_H = R * x_h
// R \in \R^{H x h}
void
CudaMgInterpolatorNested::restrict_zero(CudaVectorInterface& x_H,
                                        const CudaVectorInterface& x_h) const
{
  GlobalTimer.count("---> restrict device");
  check_matrices(x_H.n_comp, x_H.n, x_h.n);
  // x_H.zero();
  restrict_device->vmult(x_H, x_h, 1, 0);
}

// x_H -> x_h
// x_H < x_h
// x_h = R * x_H
// R \in \R^{h x H}
void
CudaMgInterpolatorNested::prolongate_add(GlobalVector& x_h,
                                         const GlobalVector& x_H) const
{
  GlobalTimer.count("---> prolongate host");
  TimePattern pattern(x_H.ncomp());
  pattern.identity();
  check_matrices(x_H.ncomp(), x_H.n(), x_h.n());
  prolongate->vmult_time(x_h, x_H, pattern, 1);
}

// x_H -> x_h
// x_H < x_h
// x_h = R * x_H
// R \in \R^{h x H}
void
CudaMgInterpolatorNested::prolongate_add(CudaVectorInterface& x_h,
                                         const CudaVectorInterface& x_H) const
{
  GlobalTimer.count("---> prolongate device");
  check_matrices(x_H.n_comp, x_H.n, x_h.n);
  prolongate_device->vmult(x_h, x_H, 1, 0);
}

}