#include "cudasolver.h"

#include <assert.h>
#include <iostream>
#include <stdlib.h>

#include "check_cuda.h"
#include "cublas_v2.h"
#include "cudacsrmatrixinterface.h"
#include "cusparse.h"
#include "gascoigne.h"
#include "matrixinterface.h"
#include "pressurefilter.h"
#include "solverdata.h"
#include "sparseblockmatrix.h"
#include "stopwatch.h"

#include "columnstencil.h"

namespace Gascoigne {
extern Timer GlobalTimer;
template<int N>
class FMatrixBlock;

CudaSolver::CudaSolver()
{
  CHECK_CUSPARSE(cusparseCreate(&sparse_handle));
  CHECK_CUBLAS(cublasCreate(&blas_handle));

  hn_zero = std::make_shared<SimpleMatrix>();
  hn_average = std::make_shared<SimpleMatrix>();
  hn_distribute = std::make_shared<SimpleMatrix>();
}

CudaSolver::~CudaSolver()
{
  cusparseDestroy(sparse_handle);
}

void
CudaSolver::ActivateCuda(std::initializer_list<const Vector*> vectors) const
{
  GlobalTimer.count("---> ActivateCuda");
  GlobalTimer.start("---> copy");
  for (auto vec : vectors) {
    InitCV(*vec);
  }
  GlobalTimer.stop("---> copy");
  GlobalTimer.start("--> oncuda");
  oncuda = true;
}

void
CudaSolver::DeactivateCuda(std::initializer_list<Vector*> vectors) const
{
  GlobalTimer.count("---> DeactivateCuda");
  GlobalTimer.stop("--> oncuda");
  GlobalTimer.start("---> copy");
  for (auto vec : vectors) {
    CopyBack(*vec);
  }
  GlobalTimer.stop("---> copy");
  oncuda = false;
}

CudaVectorInterface&
CudaSolver::InitCV(const Vector& u) const
{
  GlobalTimer.count("---> InitCV");
  if (cva.find(u) == cva.end()) {
    cva.emplace(u, new CudaVectorInterface(GetGV(u)));
  } else {
    *cva[u] = GetGV(u);
  }
  return *cva[u];
}

CudaVectorInterface&
CudaSolver::GetCV(const Vector& u) const
{
  return *cva[u];
}

GlobalVector&
CudaSolver::CopyBack(Vector& gu) const
{
  GlobalTimer.count("---> CopyBack");
  GlobalVector& u = GetGV(gu);
  cva[gu]->copy_back(u);
  return u;
}

void
CudaSolver::DeleteVector(Vector& p) const
{
  StdSolver::DeleteVector(p);
  cva.Delete(p);
}

void
CudaSolver::BasicInit(const ParamFile& paramfile, const int dimension)
{
  StdSolver::BasicInit(paramfile, dimension);
}

MatrixInterface*
CudaSolver::NewMatrix(int ncomp, const std::string& matrixtype)
{
  this->ncomp = ncomp;
  this->matrixtype = matrixtype;
  return StdSolver::NewMatrix(ncomp, matrixtype);
}

void
CudaSolver::AssembleMatrix(Matrix& A, const Vector& u, double d) const
{
  StdSolver::AssembleMatrix(A, u, d);
  GlobalTimer.start("---> matrix");
  const MatrixInterface& sbm = GetMatrix(A);

  if ((matrixtype == "block") || (matrixtype == "sparseumf") ||
      (matrixtype == "vanka")) {
    mat = std::make_shared<CudaCSRMatrixInterface>(sparse_handle, sbm);
    dia_invers =
      CudaCSRMatrixInterface::get_inverse_diagonal(sparse_handle, sbm);
  }
  GlobalTimer.stop("---> matrix");
}

template<typename S>
void
build_structure(SparseStructure& saHNAverage,
                SparseStructure& saHNDistribute,
                S& hnq2,
                size_t elements)
{
  for (const auto& entry : hnq2) {
    for (size_t i = 0; i < elements; ++i) {
      saHNAverage.build_add(entry.first, entry.second[i]);
      saHNDistribute.build_add(entry.second[i], entry.first);
    }
  }
}

template<typename S>
void
set_values(std::shared_ptr<SimpleMatrix> hn_average,
           std::shared_ptr<SimpleMatrix> hn_distribute,
           S& hnq2,
           size_t elements,
           MatrixEntryType values[])
{
  for (const auto& entry : hnq2) {
    for (size_t i = 0; i < elements; ++i) {
      hn_average->GetValue(entry.first, entry.second[i]) = values[i];
      hn_distribute->GetValue(entry.second[i], entry.first) = values[i];
    }
  }
}

void
CudaSolver::NewMesh(const GascoigneMesh* mp)
{
  StdSolver::NewMesh(mp);

  nnodes = mp->nnodes();

  DiscretizationInterface* discretisation = this->GetDiscretization();
  auto disc_name = discretisation->GetName();

  SparseStructure diagonalStructure;
  diagonalStructure.build_begin(nnodes);
  for (int i = 0; i < nnodes; ++i) {
    diagonalStructure.build_add(i, i);
  }
  diagonalStructure.build_end();

  const auto& hnstructure = mp->GetHangingIndexHandler();
  const std::map<int, std::array<int, 3>>* hnq2 = hnstructure.GetStructure();
  const std::map<int, std::array<int, 9>>* hnq2face =
    hnstructure.GetStructureFace();

  SparseStructure saHNAverage;
  saHNAverage.build_begin(nnodes);
  SparseStructure saHNDistribute;
  saHNDistribute.build_begin(nnodes);

  if (disc_name == "Q12d") {
    build_structure(saHNAverage, saHNDistribute, *hnq2, 2);
  } else if (disc_name == "Q22d") {
    build_structure(saHNAverage, saHNDistribute, *hnq2, 3);
  } else if (disc_name == "Q13d") {
    build_structure(saHNAverage, saHNDistribute, *hnq2, 2);
    build_structure(saHNAverage, saHNDistribute, *hnq2face, 4);
  } else if (disc_name == "Q23d") {
    build_structure(saHNAverage, saHNDistribute, *hnq2, 3);
    build_structure(saHNAverage, saHNDistribute, *hnq2face, 9);
  }

  for (int i = 0; i < nnodes; ++i) {
    saHNDistribute.build_add(i, i);
    saHNAverage.build_add(i, i);
  }
  saHNAverage.build_end();
  saHNDistribute.build_end();

  hn_zero->ReInit(&diagonalStructure);
  hn_average->ReInit(&saHNAverage);
  hn_distribute->ReInit(&saHNDistribute);
  for (int i = 0; i < nnodes; ++i) {
    hn_zero->GetValue(i, i) = 1;
    hn_average->GetValue(i, i) = 1;
    hn_distribute->GetValue(i, i) = 1;
  }

  if (disc_name == "Q12d") {
    MatrixEntryType wei2d[] = { 0.5, 0.5 };
    set_values(hn_average, hn_distribute, *hnq2, 2, wei2d);
  } else if (disc_name == "Q22d") {
    MatrixEntryType wei2d[] = { 0.375, 0.75, -0.125 };
    set_values(hn_average, hn_distribute, *hnq2, 3, wei2d);
  } else if (disc_name == "Q13d") {
    MatrixEntryType wei2d[] = { 0.5, 0.5 };
    MatrixEntryType wei3d[] = { 0.25, 0.25, 0.25, 0.25 };
    set_values(hn_average, hn_distribute, *hnq2, 2, wei2d);
    set_values(hn_average, hn_distribute, *hnq2face, 4, wei3d);
  } else if (disc_name == "Q23d") {
    MatrixEntryType wei2d[] = { 0.375, 0.75, -0.125 };
    MatrixEntryType wei3d[9];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        wei3d[3 * i + j] = wei2d[i] * wei2d[j];
      }
    }
    set_values(hn_average, hn_distribute, *hnq2, 3, wei2d);
    set_values(hn_average, hn_distribute, *hnq2face, 9, wei3d);
  }

  for (const auto entry : *hnq2) {
    hn_distribute->GetValue(entry.first, entry.first) = 0;
    hn_zero->GetValue(entry.first, entry.first) = 0;
  }
  for (const auto entry : *hnq2face) {
    hn_distribute->GetValue(entry.first, entry.first) = 0;
    hn_zero->GetValue(entry.first, entry.first) = 0;
  }
}

void
CudaSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  StdSolver::SetProblem(PDX);

  ncomp = PDX.GetNcomp();

  SparseStructure diagonalStructure;
  diagonalStructure.build_begin(nnodes);
  for (int i = 0; i < nnodes; ++i) {
    diagonalStructure.build_add(i, i);
  }
  diagonalStructure.build_end();

  TimePattern identity(ncomp);
  identity.identity();
  hn_distribute_device = std::make_shared<CudaCSRMatrixInterface>(
    sparse_handle, hn_distribute, &identity, 0);
  hn_zero_device = std::make_shared<CudaCSRMatrixInterface>(
    sparse_handle, hn_zero, &identity, 0);
  hn_average_device = std::make_shared<CudaCSRMatrixInterface>(
    sparse_handle, hn_average, &identity, 0);

  const BoundaryIndexHandler& BIH = GetMesh()->GetBoundaryIndexHandler();
  const DirichletData* DD = GetProblemDescriptor()->GetDirichletData();

  ColumnStencil CS;
  CS.memory(&diagonalStructure);

  std::vector<TimePattern> values(nnodes, identity);
  for (int col : DD->dirichlet_colors()) {

    nvector<int> comps = DD->components_on_color(col);

    TimePattern pattern(ncomp);
    pattern.identity();
    for (int i : comps) {
      pattern(i, i) = 0;
    }

    nvector<int> nodes = BIH.Verteces(col);
    for (int i : nodes) {
      for (int j = 0; j < ncomp; ++j) {
        values[i](j, j) *= pattern(j, j);
      }
    }
  }
  dirichlet_zeros = std::make_shared<CudaCSRMatrixInterface>(
    sparse_handle, values, &CS, nnodes);
}

void
CudaSolver::vmult(const Matrix& A, Vector& gy, const Vector& gx, double d) const
{
  if (!oncuda) {
    GlobalTimer.start("---> vmult");
    ActivateCuda({ &gy, &gx });
    GlobalTimer.stop("---> vmult");
    vmult(A, gy, gx, d);
    GlobalTimer.start("---> vmult");
    DeactivateCuda({ &gy });
    GlobalTimer.stop("---> vmult");
    // StdSolver::vmult(A, gy, gx, d);
    return;
  }
  GlobalTimer.start("---> vmult");
  CudaVectorInterface& y = GetCV(gy);
  CudaVectorInterface& x = GetCV(gx);
  mat->vmult(y, x, d, 1);
  GlobalTimer.stop("---> vmult");
}

void
CudaSolver::vmulteq(const Matrix& A,
                    Vector& gy,
                    const Vector& gx,
                    double d) const
{
  if (!oncuda) {
    GlobalTimer.start("---> vmult");
    ActivateCuda({ &gy, &gx });
    GlobalTimer.stop("---> vmult");
    vmulteq(A, gy, gx, d);
    GlobalTimer.start("---> vmult");
    DeactivateCuda({ &gy });
    GlobalTimer.stop("---> vmult");
    // StdSolver::vmulteq(A, gy, gx, d);
    return;
  }
  GlobalTimer.start("---> vmult");
  Zero(gy);
  GlobalTimer.stop("---> vmult");
  vmult(A, gy, gx, d);
}

void
CudaSolver::smooth(int niter,
                   const Matrix& A,
                   Vector& x,
                   const Vector& y,
                   Vector& h) const
{
  if (!oncuda) {
    StdSolver::smooth(niter, A, x, y, h);
    return;
  }
  GlobalTimer.start("---> smooth");
  double omega = GetSolverData().GetOmega();
  // std::cout << Norm(h) << " " << Norm(x) << std::endl;
  for (int iter = 0; iter < niter; iter++) {
    if (GetSolverData().GetLinearSmooth() == "jacobi") {
      GlobalTimer.stop("---> smooth");
      MatrixResidual(A, h, x, y);
      GlobalTimer.start("---> smooth");
      CudaVectorInterface tmp(GetCV(h).n, GetCV(h).n_comp);
      dia_invers->vmult(tmp, GetCV(h), 1, 0);
      GetCV(h) = tmp;
      Add(x, omega, h);
    } else if (GetSolverData().GetLinearSmooth() == "richardson") {
      MatrixResidual(A, h, x, y);
      Add(x, omega, h);
    }
    // SubtractMean(gx); // TODO: Implement for cuda
  }
  GlobalTimer.stop("---> smooth");
}

/**
 * gy = gb - A gx
 */
void
CudaSolver::MatrixResidual(const Matrix& A,
                           Vector& gy,
                           const Vector& gx,
                           const Vector& gb) const
{
  if (!oncuda) {
    StdSolver::MatrixResidual(A, gy, gx, gb);
    return;
  }

  Equ(gy, 1, gb);
  vmult(A, gy, gx, -1); // y += -1 * A * x

  if (GetPfilter().Active()) {
    DeactivateCuda({ &gy });
    SubtractMeanAlgebraic(gy);
    ActivateCuda({ &gy });
  }
}

/**
 * dst = 0
 */
void
CudaSolver::Zero(Vector& dst) const
{
  if (!oncuda) {
    StdSolver::Zero(dst);
    return;
  }
  GetCV(dst).zero();
}

/**
 * gdst = s * gsrc
 */
void
CudaSolver::Equ(Vector& gdst, double s, const Vector& gsrc) const
{
  if (!oncuda) {
    StdSolver::Equ(gdst, s, gsrc);
    return;
  }
  CudaVectorInterface& dst = GetCV(gdst);
  dst = GetCV(gsrc);
  if (s != 1) {
    dst.scal(blas_handle, s);
  }
}

void
CudaSolver::Add(Vector& gdst, double s, const Vector& gsrc) const
{
  if (!oncuda) {
    StdSolver::Add(gdst, s, gsrc);
    return;
  }
  GetCV(gdst).add(blas_handle, s, GetCV(gsrc));
}

double
CudaSolver::Norm(const Vector& gdst) const
{
  if (!oncuda) {
    return StdSolver::Norm(gdst);
  }
  return GetCV(gdst).norm(blas_handle);
}

void
CudaSolver::HNZero(const Vector& x) const
{
  if (!IsOnCuda()) {
    StdSolver::HNZero(x);
    return;
  }

  CudaVectorInterface tmp(GetCV(x).n, GetCV(x).n_comp);
  hn_zero_device->vmult(tmp, GetCV(x), 1, 0);
  GetCV(x) = tmp;
}

void
CudaSolver::HNDistribute(Vector& x) const
{
  if (!IsOnCuda()) {
    StdSolver::HNDistribute(x);
    return;
  }
  CudaVectorInterface& x_cuda = GetCV(x);

  CudaVectorInterface tmp(x_cuda.n, x_cuda.n_comp);
  // tmp = x_cuda;
  hn_distribute_device->vmult(tmp, x_cuda, 1, 0);
  x_cuda = tmp;
}

void
CudaSolver::HNAverage(const Vector& x) const
{
  if (!IsOnCuda()) {
    StdSolver::HNAverage(x);
    return;
  }

  CudaVectorInterface tmp(GetCV(x).n, GetCV(x).n_comp);
  // tmp = GetCV(x);
  hn_average_device->vmult(tmp, GetCV(x), 1, 0);
  GetCV(x) = tmp;
}

void
CudaSolver::SetBoundaryVectorZero(Vector& gf) const
{
  if (!IsOnCuda()) {
    StdSolver::SetBoundaryVectorZero(gf);
    return;
  }
  CudaVectorInterface tmp(GetCV(gf).n, GetCV(gf).n_comp);
  dirichlet_zeros->vmult(tmp, GetCV(gf), 1, 0);
  GetCV(gf) = tmp;
}

} // namespace Gascoigne
