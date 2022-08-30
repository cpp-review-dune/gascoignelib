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
ActivateCuda(CudaMultiLevelSolver* smls,
             IndexType finelevel,
             IndexType coarselevel,
             std::initializer_list<const Vector*> vectors)
{
  for (size_t i = coarselevel; i <= finelevel; i++) {
    static_cast<CudaSolver*>(smls->GetSolver(static_cast<IndexType>(i)))
      ->ActivateCuda(vectors);
  }
}

void
DeactivateCuda(CudaMultiLevelSolver* smls,
               IndexType finelevel,
               IndexType coarselevel,
               std::initializer_list<Vector*> vectors)
{
  for (size_t i = coarselevel; i <= finelevel; i++) {
    static_cast<CudaSolver*>(smls->GetSolver(static_cast<IndexType>(i)))
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

auto CudaMultiLevelSolver::NewSolver(IndexType /*solverlevel*/) -> StdSolver*
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

  for (IndexType l = 0; l < nlevels() - 1; ++l) {
    _Interpolator[l] = new CudaMgInterpolatorNested;
  }

  //
  // Interpolator [l] :   interpoliert   (l+1)->l  (fein->grob)
  //
  for (IndexType level = 0; level < nlevels() - 1; ++level) {
    IndexType sl = nlevels() - level - 2;

    const MeshTransferInterface* MT = GetMeshAgent()->GetTransfer(sl);
    assert(MT);
    assert(_Interpolator[level]);
    GetSolver(level)->ConstructInterpolator(_Interpolator[level], MT);
  }
}

void
CudaMultiLevelSolver::LinearMg(IndexType finelevel,
                               IndexType coarselevel,
                               const Matrix& A,
                               Vector& u,
                               const Vector& f,
                               CGInfo& info)
{
  GlobalTimer.start("--> LinearMg");

  for (IndexType level = coarselevel; level < nlevels(); level++) {
    if (GetSolver(level)->DirectSolver()) {
      coarselevel = level;
    }
  }
  // wir haben auf einem hoeheren level einen direkten loeser...
  if (coarselevel > finelevel) {
    coarselevel = finelevel;
  }

  assert(finelevel >= coarselevel);

  IndexType nl = nlevels();
  DoubleVector res(nl, 0.), rw(nl, 0.);

  Vector mg0("mg0_");
  Vector mg1("mg1_");
  ReInitVector(mg0);
  ReInitVector(mg1);

  ActivateCuda(this, finelevel, coarselevel, { &mg0, &f, &mg1, &u });

  CudaSolver* solver = static_cast<CudaSolver*>(GetSolver(finelevel));

  solver->Equ(mg0, 1., f);

  solver->MatrixResidual(A, mg1, u, mg0);

  res[finelevel] = solver->Norm(mg1);

  rw[finelevel] = 0.;
  info.check(res[finelevel], rw[finelevel]);

  bool reached = false; // mindestens einen schritt machen
  for (IndexType it = 0; !reached; it++) {
    std::string p = DataP.MgType();
    std::string p0 = p;
    if (p == "F") {
      p0 = "W";
    }
    mgstep(res, rw, finelevel, finelevel, coarselevel, p0, p, A, u, mg0, mg1);
    reached = info.check(res[finelevel], rw[finelevel]);
  }

  DeactivateCuda(this, finelevel, coarselevel, { &u });

  DeleteVector(mg0);
  DeleteVector(mg1);
  GlobalTimer.stop("--> LinearMg");
}

/**
 * Expects vor GetSolver(l) to be Cuda Activated
 */
void
CudaMultiLevelSolver::mgstep(std::vector<double>& res,
                             std::vector<double>& rw,
                             IndexType l,
                             IndexType finelevel,
                             IndexType coarselevel,
                             std::string& p0,
                             std::string p,
                             const Matrix& A,
                             Vector& u,
                             Vector& b,
                             Vector& v)
{
  GlobalTimer.count("--> mgstep");

  CudaSolver* solver_l = static_cast<CudaSolver*>(GetSolver(l));
  CudaSolver* solver_lm1 = static_cast<CudaSolver*>(GetSolver(l - 1));

  if (l == coarselevel) {
    if (p == "F") {
      p0 = "V";
    }
    solver_l->smooth_exact(A, u, b, v);
    if (coarselevel == finelevel) {
      solver_l->MatrixResidual(A, v, u, b);
      res[l] = solver_l->Norm(v);
    }
    return;
  }

  solver_l->smooth_pre(A, u, b, v);
  solver_l->MatrixResidual(A, v, u, b);

  dynamic_cast<CudaMgInterpolatorNested*>(_Interpolator[l - 1])
    ->restrict_zero(solver_lm1->GetCV(b), solver_l->GetCV(v));

  solver_lm1->HNDistribute(b);
  solver_lm1->SetBoundaryVectorZero(b);

  solver_lm1->Zero(u);

  IndexType j = 0;
  if (p0 == "V") {
    j = 1;
  } else if (p0 == "W") {
    j = 2;
  } else if (p0 == "F") {
    j = 3;
  }

  for (IndexType i = 0; i < j; ++i) {
    mgstep(res, rw, l - 1, finelevel, coarselevel, p0, p, A, u, b, v);
  }

  if ((l == 0) && (p == "F")) {
    p0 = "W";
  }

  rw[l] = solver_lm1->Norm(u);
  solver_lm1->HNAverage(u);

  solver_l->Zero(v);

  dynamic_cast<CudaMgInterpolatorNested*>(_Interpolator[l - 1])
    ->prolongate_add(solver_l->GetCV(v), solver_lm1->GetCV(u));

  solver_lm1->HNZero(u);

  solver_l->HNZero(v);
  solver_l->SetBoundaryVectorZero(v);

  solver_l->Add(u, DataP.MgOmega(), v);

  solver_l->smooth_post(A, u, b, v);
  solver_l->MatrixResidual(A, v, u, b);
  res[l] = solver_l->Norm(v);
}

/*-------------------------------------------------------------*/
} // namespace Gascoigne
