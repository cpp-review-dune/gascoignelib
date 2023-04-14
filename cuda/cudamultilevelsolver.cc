#include "cudamultilevelsolver.h"

#include <assert.h>
#include <string>
#include <vector>

#include "cginfo.h"
#include "cudamgintepolatornested.h"
#include "cudasolver.h"
#include "gascoigne.h"
#include "mginterpolatorinterface.h"
#include "stdmultilevelsolverdata.h"
#include "stdsolver.h"
#include "stopwatch.h"

#include "cg.h"
#include "compose_name.h"
#include "gascoignemultigridmesh.h"
#include "gmres.h"
#include "mginterpolatormatrix.h"
#include "mginterpolatornested.h"
#include "stdmultilevelsolver.h"
#include <iomanip>

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne {
extern Timer GlobalTimer;
// **************************************************

void
CudaMultiLevelSolver::ActivateCuda(std::initializer_list<const Vector*> vectors)
{
  for (size_t i = 0; i <= FinestLevel(); i++) {
    static_cast<CudaSolver*>(this->GetSolver(static_cast<int>(i)))
      ->ActivateCuda(vectors);
  }
}

void
CudaMultiLevelSolver::DeactivateCuda(std::initializer_list<Vector*> vectors)
{
  for (size_t i = 0; i <= FinestLevel(); i++) {
    static_cast<CudaSolver*>(this->GetSolver(static_cast<int>(i)))
      ->DeactivateCuda(vectors);
  }
}

CudaMultiLevelSolver::CudaMultiLevelSolver()
  : StdMultiLevelSolver()
{}

CudaMultiLevelSolver::CudaMultiLevelSolver(const MeshAgent* MAP,
                                           const ParamFile& paramfile,
                                           const ProblemContainer* PC,
                                           const FunctionalContainer* FC)
  : StdMultiLevelSolver(MAP, paramfile, PC, FC)
{}

CudaMultiLevelSolver::~CudaMultiLevelSolver() {}

auto
CudaMultiLevelSolver::NewSolver(int /*solverlevel*/) -> StdSolver*
{
  return new CudaSolver;
}

void
CudaMultiLevelSolver::NewMgInterpolator()
{
  if (DataP.LinearSolve() == "direct")
    return; // No interpolator for direct solver

  for (size_t i = 0; i < _Interpolator.size(); ++i) {
    assert(_Interpolator[i] != NULL);
    delete _Interpolator[i];
    _Interpolator[i] = NULL;
  }
  _Interpolator.resize(nlevels() - 1, NULL);

  for (int l = 0; l < nlevels() - 1; ++l) {
    _Interpolator[l] = new CudaMgInterpolatorNested;
  }

  //
  // Interpolator [l] :   interpoliert   (l+1)->l  (fein->grob)
  //
  for (int level = 0; level < nlevels() - 1; ++level) {
    int sl = nlevels() - level - 2;

    const MeshTransferInterface* MT = GetMeshAgent()->GetTransfer(sl);
    assert(MT);
    assert(_Interpolator[level]);
    GetSolver(level)->ConstructInterpolator(_Interpolator[level], MT);
  }
}

void
CudaMultiLevelSolver::NewtonLinearSolve(const Matrix& A,
                                        Vector& x,
                                        const Vector& b,
                                        CGInfo& info)
{
  _clock_solve.start();
  info.reset();
  GetSolver(FinestLevel())->Zero(x);

  assert(DataP.LinearSolve() == "mg" || DataP.LinearSolve() == "gmres" ||
         DataP.LinearSolve() == "direct");

  if (DataP.LinearSolve() == "mg") {
    int clevel = std::max(DataP.CoarseLevel(), 0);
    if (DataP.CoarseLevel() == -1)
      clevel = FinestLevel();
    ActivateCuda({ &x, &b });
    LinearMg(ComputeLevel, clevel, A, x, b, info);
    DeactivateCuda({ &x });
  } else if (DataP.LinearSolve() == "direct") {
    int clevel = FinestLevel();
    ActivateCuda({ &x, &b });
    LinearMg(ComputeLevel, clevel, A, x, b, info);
    DeactivateCuda({ &x });
  } else if (DataP.LinearSolve() == "gmres") {
    ActivateCuda({ &x, &b });
    Gmres(A, x, b, info);
    DeactivateCuda({ &x });
  }

  GetSolver(ComputeLevel)->SubtractMean(x);
  _clock_solve.stop();
}

/**
 * Multigrid interpolation of vector v to coarse level vector b
 *
 * @param level of outputvector b
 */
void
CudaMultiLevelSolver::RestrictZero(IndexType level,
                                   Vector& b,
                                   const Vector& v) const
{
  dynamic_cast<CudaMgInterpolatorNested*>(_Interpolator[level])
    ->restrict_zero(
      static_cast<const CudaSolver*>(GetSolver(level))->GetCV(b),
      static_cast<const CudaSolver*>(GetSolver(level + 1))->GetCV(v));
}

/**
 * Multigrid interpolation of coarse level vector v to fine level vector b
 *
 * @param level of outputvector b
 */
void
CudaMultiLevelSolver::ProlongateAdd(IndexType level,
                                    Vector& b,
                                    const Vector& v) const
{
  dynamic_cast<CudaMgInterpolatorNested*>(_Interpolator[level])
    ->prolongate_add(
      static_cast<const CudaSolver*>(GetSolver(level + 1))->GetCV(b),
      static_cast<const CudaSolver*>(GetSolver(level))->GetCV(v));
}

/*-------------------------------------------------------------*/
} // namespace Gascoigne
