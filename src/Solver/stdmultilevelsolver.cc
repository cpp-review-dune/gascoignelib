/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011 by the Gascoigne 3D
 *authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/

#include "stdmultilevelsolver.h"
#include "cg.h"
#include "compose_name.h"
#include "gascoignemultigridmesh.h"
#include "gmres.h"
#include "mginterpolatormatrix.h"
#include "mginterpolatornested.h"
#include <iomanip>

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne {
extern Timer GlobalTimer;

StdMultiLevelSolver::~StdMultiLevelSolver()
{
  // ViewProtocoll();

  for (int i = 0; i < GetSolverPointers().size(); i++) {
    if (GetSolverPointers()[i]) {
      delete GetSolverPointers()[i];
      GetSolverPointers()[i] = NULL;
    }
  }
  for (int i = 0; i < _Interpolator.size(); i++) {
    if (_Interpolator[i]) {
      delete _Interpolator[i];
      _Interpolator[i] = NULL;
    }
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::ViewProtocoll() const
{
  // old...
  GlobalTimer.print100();
  assert(0);
  // cout <<
  // "\n************************************************************************\n\n";
  // cout << "StdMultiLevelSolver\tTIME\n";
  // cout << "  Residual\t\t" << _clock_residual.read() << endl;
  // cout << "  Solve\t\t\t" << _clock_solve.read() << endl;
  // cout << "MATRIX\t\t\tTIME\n";

  // double vm=0.;
  // double il=0.;
  // double so=0.;
  // double ca=0., ci=0., cs=0.;
  // double re=0;
  // for(int level=0;level<GetSolverPointers().size();level++)
  //   {
  // 	const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(level));
  // 	assert(S);
  // 	vm += S->clock_vmult();
  // 	il += S->clock_ilu();
  // 	so += S->clock_solve();
  // 	ca += S->clock_computematrix();
  // 	ci += S->clock_computeilu();
  // 	cs += S->clock_computesolver();
  // 	re += S->clock_residual();
  //   }

  // cout << "  vmult\t\t\t" << vm << endl;
  // cout << "  smooth\t\t\t" << il << endl;
  // cout << "  solve\t\t\t" << so << endl;
  // cout << "  compute matrix\t" << ca << endl;
  // cout << "  compute ilu\t\t" << ci << endl;
  // cout << "  compute solver\t" << cs << endl;
  // cout << "VECTOR\t\t\tTIME\n";
  // cout << "  residual\t\t" << re << endl;
  // cout <<
  // "\n************************************************************************\n";
}

/*-------------------------------------------------------------*/

StdMultiLevelSolver::StdMultiLevelSolver()
  : _MAP(NULL)
  , oldnlevels(-1)
  , MON()
  , _PD(0)
  , _PC(0)
  , _FC(0)
{
}

StdMultiLevelSolver::StdMultiLevelSolver(const MeshAgent* MAP,
                                         const ParamFile& paramfile,
                                         const ProblemContainer* PC,
                                         const FunctionalContainer* FC)
  : _MAP(MAP)
  , oldnlevels(-1)
  , _paramfile(paramfile)
  , MON(_paramfile, 1)
  , _PD(0)
  , _PC(PC)
  , _FC(FC)
  , DataP(_paramfile)
{
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::BasicInit(const MeshAgent* MAP,
                               const ParamFile& paramfile,
                               const ProblemContainer* PC,
                               const FunctionalContainer* FC)
{
  _MAP = MAP;

  _paramfile = paramfile;

  SetProblemContainer(PC);
  SetFunctionalContainer(FC);
  DataP = StdMultiLevelSolverData();
  DataP.BasicInit(_paramfile);
  MON = Monitor(_paramfile, 1);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::SetProblem(const std::string& label)
{
  _PD = GetProblemContainer()->GetProblem(label);
  for (int level = 0; level < nlevels(); ++level) {
    int solverlevel = nlevels() - 1 - level;
    assert(GetSolver(solverlevel));
    GetSolver(solverlevel)->SetProblem(*_PD);
  }
}

/*-------------------------------------------------------------*/

StdSolver*
StdMultiLevelSolver::NewSolver(int solverlevel)
{
  if (DataP.Solver() == "instat") {
    std::cout << "StdTimeSolver muss an neue Matrix angepasst werden"
              << std::endl;
    abort();
    //    return new StdTimeSolver;
  } else {
    return new StdSolver;
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::ReInitMatrix(const Matrix& A)
{
  for (int level = 0; level < nlevels(); ++level) {
    GetSolver(level)->ReInitMatrix(A);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::ReInitVector(Vector& v)
{
  for (int level = 0; level < nlevels(); ++level) {
    GetSolver(level)->ReInitVector(v);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::ReInitVector(Vector& v, int comp)
{
  for (int level = 0; level < nlevels(); ++level) {
    GetSolver(level)->ReInitVector(v, comp);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::NewSolvers()
{
  oldnlevels = GetSolverPointers().size();

  if (oldnlevels > nlevels()) {
    for (int l = oldnlevels - 1; l >= nlevels(); l--) {
      delete GetSolverPointers()[l];
      GetSolverPointers()[l] = NULL;
    }
  }
  GetSolverPointers().resize(nlevels(), NULL);
  ComputeLevel = GetSolverPointers().size() - 1;

  for (int level = 0; level < nlevels(); ++level) {
    int solverlevel = nlevels() - 1 - level;

    // new Solvers
    if (GetSolver(solverlevel) == NULL) {
      GetSolverPointer(solverlevel) = NewSolver(solverlevel);
      GetSolver(solverlevel)
        ->BasicInit(_paramfile, GetMeshAgent()->GetDimension());
    }
  }
}

/*-------------------------------------------------------------*/

// void StdMultiLevelSolver::RegisterMatrix()
// {
//   for (int level = 0; level < nlevels(); ++level)
//   {
//     GetSolver(level)->RegisterMatrix();
//   }
// }

// void StdMultiLevelSolver::ReInitMatrix()
// {
//   for (int level = 0; level < nlevels(); ++level)
//   {
//     GetSolver(level)->ReInitMatrix();
//   }
// }

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::SolverNewMesh()
{
  for (int level = 0; level < nlevels(); ++level) {
    const GascoigneMesh* MIP = GetMeshAgent()->GetMesh(level);
    assert(MIP);

    int solverlevel = nlevels() - 1 - level;
    GetSolver(solverlevel)->NewMesh(MIP);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::NewMgInterpolator()
{
  if (DataP.LinearSolve() == "direct")
    return; // No interpolator for direct solver

  for (int i = 0; i < _Interpolator.size(); ++i) {
    assert(_Interpolator[i] != NULL);
    delete _Interpolator[i];
    _Interpolator[i] = NULL;
  }
  _Interpolator.resize(nlevels() - 1, NULL);

  for (int l = 0; l < nlevels() - 1; ++l)
    //  _Interpolator[l] = new MgInterpolatorMatrix;
    _Interpolator[l] = new MgInterpolatorNested;

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

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::ReInit()
{
  DataP.CountResidual() = 0;
  NewSolvers();
  SolverNewMesh();
  NewMgInterpolator();
}

/*-------------------------------------------------------------*/

const std::vector<std::string>
StdMultiLevelSolver::GetFunctionalNames() const
{
  std::vector<std::string> names;
  if (!GetFunctionalContainer())
    return names;

  for (FunctionalContainer::const_iterator it =
         GetFunctionalContainer()->begin();
       it != GetFunctionalContainer()->end();
       ++it)
    names.push_back(it->second->GetName());
  return names;
}

/*-------------------------------------------------------------*/

const DoubleVector
StdMultiLevelSolver::GetExactValues() const
{
  if (!GetFunctionalContainer())
    return DoubleVector(0);

  int n = GetFunctionalContainer()->size();
  DoubleVector j(n, 0.);
  int i = 0;
  for (FunctionalContainer::const_iterator it =
         GetFunctionalContainer()->begin();
       it != GetFunctionalContainer()->end();
       ++it, ++i)
    j[i] = it->second->ExactValue();
  return j;
}

/*-------------------------------------------------------------*/

const DoubleVector
StdMultiLevelSolver::ComputeFunctionals(Vector& f,
                                        const Vector& u,
                                        FunctionalContainer* FC)
{
  if (!FC)
    return DoubleVector(0);
  int n = FC->size();
  DoubleVector j(n, 0.);
  int i = 0;
  for (FunctionalContainer::const_iterator it = FC->begin(); it != FC->end();
       ++it, ++i) {
    j[i] = GetSolver(ComputeLevel)->ComputeFunctional(f, u, it->second);
  }
  return j;
}

/*-------------------------------------------------------------*/

const DoubleVector
StdMultiLevelSolver::ComputeFunctionals(Vector& f, const Vector& u)
{
  if (!GetFunctionalContainer())
    return DoubleVector(0);
  int n = GetFunctionalContainer()->size();
  DoubleVector j(n, 0.);
  int i = 0;
  for (FunctionalContainer::const_iterator it =
         GetFunctionalContainer()->begin();
       it != GetFunctionalContainer()->end();
       ++it, ++i) {
    j[i] = GetSolver(ComputeLevel)->ComputeFunctional(f, u, it->second);
  }
  return j;
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::InterpolateSolution(Vector& u,
                                         const GlobalVector& uold) const
{
  if (oldnlevels <= 0)
    return;

  GetSolver(FinestLevel())->InterpolateSolution(u, uold);

  for (int l = 0; l < FinestLevel(); l++) {
    GetSolver(l)->Zero(u);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::vmulteq(const Matrix& A, Vector& y, const Vector& x) const
{
  GetSolver(ComputeLevel)->vmulteq(A, y, x, 1.);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::LinearMg(int finelevel,
                              int coarselevel,
                              const Matrix& A,
                              Vector& u,
                              const Vector& f,
                              CGInfo& info)
{
  int clevel = coarselevel;
  for (int level = coarselevel; level < nlevels(); level++) {
    if (GetSolver(level)->DirectSolver())
      clevel = level;
  }
  // wir haben auf einem hoeheren level einen direkten loeser...
  if (clevel > finelevel) {
    clevel = finelevel;
  }

  assert(finelevel >= clevel);

  int nl = nlevels();
  DoubleVector res(nl, 0.), rw(nl, 0.);

  Vector mg0("_mg0_");
  Vector mg1("_mg1_");
  ReInitVector(mg0);
  ReInitVector(mg1);

  GetSolver(finelevel)->Equ(mg0, 1., f);

  GetSolver(finelevel)->MatrixResidual(A, mg1, u, mg0);
  res[finelevel] = GetSolver(finelevel)->Norm(mg1);
  rw[finelevel] = 0.;
  info.check(res[finelevel], rw[finelevel]);

  bool reached = false; // mindestens einen schritt machen
  for (int it = 0; !reached; it++) {
    string p = DataP.MgType();
    string p0 = p;
    if (p == "F")
      p0 = "W";
    mgstep(res, rw, finelevel, finelevel, clevel, p0, p, A, u, mg0, mg1);
    reached = info.check(res[finelevel], rw[finelevel]);
  }

  DeleteVector(mg0);
  DeleteVector(mg1);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::mgstep(vector<double>& res,
                            vector<double>& rw,
                            int l,
                            int finelevel,
                            int coarselevel,
                            string& p0,
                            string p,
                            const Matrix& A,
                            Vector& u,
                            Vector& b,
                            Vector& v)
{
  GlobalTimer.count("--> mgstep");

  if (l == coarselevel) {
    if (p == "F") {
      p0 = "V";
    }
    GetSolver(l)->smooth_exact(A, u, b, v);
    if (l == finelevel) {
      GetSolver(l)->MatrixResidual(A, v, u, b);
      res[l] = GetSolver(l)->Norm(v);
    }
  } else {
    GetSolver(l)->smooth_pre(A, u, b, v);
    GetSolver(l)->MatrixResidual(A, v, u, b);

    _Interpolator[l - 1]->restrict_zero(GetSolver(l - 1)->GetGV(b),
                                        GetSolver(l)->GetGV(v));
    GetSolver(l - 1)->HNDistribute(b);
    GetSolver(l - 1)->SetBoundaryVectorZero(b);

    GetSolver(l - 1)->Zero(u);

    int j = 0;
    if (p0 == "V")
      j = 1;
    if (p0 == "W")
      j = 2;
    if (p0 == "F")
      j = 3;

    for (int i = 0; i < j; i++) {
      mgstep(res, rw, l - 1, finelevel, coarselevel, p0, p, A, u, b, v);
    }
    if ((l == 0) && (p == "F")) {
      p0 = "W";
    }
    rw[l] = GetSolver(l - 1)->Norm(u);

    GetSolver(l)->Zero(v);

    GetSolver(l - 1)->HNAverage(u);

    _Interpolator[l - 1]->prolongate_add(GetSolver(l)->GetGV(v),
                                         GetSolver(l - 1)->GetGV(u));

    GetSolver(l - 1)->HNZero(u);

    GetSolver(l)->HNZero(v);
    GetSolver(l)->SetBoundaryVectorZero(v);

    GetSolver(l)->Add(u, DataP.MgOmega(), v);

    GetSolver(l)->smooth_post(A, u, b, v);
    GetSolver(l)->MatrixResidual(A, v, u, b);
    res[l] = GetSolver(l)->Norm(v);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::Cg(Vector& x, const Vector& f, CGInfo& info)
{
  std::cerr << "\"StdMultiLevelSolver::Cg\" not written!" << std::endl;
  abort();
  //   CG<SolverInterface,StdMultiLevelSolver,GhostVector>
  //   cg(*(GetSolver(ComputeLevel)),*this);

  //   cg.solve(x,f,info);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::newton(Matrix& A,
                            Vector& u,
                            const Vector& f,
                            Vector& r,
                            Vector& w,
                            NLInfo& info)
{
  info.reset();
  double rr = NewtonResidual(r, u, f);
  bool reached = info.check(0, rr, 0.);
  NewtonOutput(info);
  NewtonPreProcess(u, f, info);
  for (int it = 1; !reached; it++) {
    NewtonMatrixControl(A, u, info);
    NewtonVectorZero(w);
    NewtonLinearSolve(A, w, r, info.GetLinearInfo());
    double rw = NewtonUpdate(rr, u, w, r, f, info);
    reached = info.check(it, rr, rw);
    NewtonOutput(info);
  }
  NewtonPostProcess(u, f, info);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::NewtonOutput(NLInfo& nlinfo) const
{
  MON.nonlinear_step(nlinfo.GetLinearInfo(), nlinfo);
}
/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::NewtonPreProcess(Vector& u,
                                      const Vector& f,
                                      NLInfo& info) const
{
  ;
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::NewtonPostProcess(Vector& u,
                                       const Vector& f,
                                       NLInfo& info) const
{
  ;
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::NewtonVectorZero(Vector& w) const
{
  GetSolver(FinestLevel())->Zero(w);
}

/*-------------------------------------------------------------*/

double
StdMultiLevelSolver::NewtonResidual(Vector& y,
                                    const Vector& x,
                                    const Vector& b) const
{
  _clock_residual.start();
  DataP.CountResidual()++;
  GetSolver(ComputeLevel)->Equ(y, 1., b);
  GetSolver(ComputeLevel)->Form(y, x, -1.);
  GetSolver(ComputeLevel)->SetPeriodicVectorZero(y);
  GetSolver(ComputeLevel)->SetBoundaryVectorZero(y);
  GetSolver(ComputeLevel)->SubtractMeanAlgebraic(y);
  _clock_residual.stop();
  return NewtonNorm(y);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::NewtonMatrixControl(Matrix& A, Vector& u, NLInfo& nlinfo)
{
  MON.new_matrix() = 0;

  int nm1 = nlinfo.control().newmatrix();
  int nm2 = nlinfo.control().matrixmustbebuild();

  if (nm1 + nm2 == 0)
    return;

  MON.new_matrix() = 1;

  AssembleMatrix(A, u, nlinfo);
  ComputeIlu(A, u);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::NewtonLinearSolve(const Matrix& A,
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
    LinearMg(ComputeLevel, clevel, A, x, b, info);
  } else if (DataP.LinearSolve() == "direct") {
    int clevel = FinestLevel();
    LinearMg(ComputeLevel, clevel, A, x, b, info);
  } else if (DataP.LinearSolve() == "gmres") {
    Gmres(A, x, b, info);
  }

  GetSolver(ComputeLevel)->SubtractMean(x);
  _clock_solve.stop();
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::Gmres(const Matrix& A,
                           Vector& x,
                           const Vector& f,
                           CGInfo& info)
{
  int n = DataP.GmresMemSize();

  StdSolver* S = dynamic_cast<StdSolver*>(GetSolver(ComputeLevel));
  GMRES gmres(*S, *this, n);

  gmres.solve(A, x, f, info);
}

/*-------------------------------------------------------------*/

double
StdMultiLevelSolver::NewtonUpdate(double& rr,
                                  Vector& x,
                                  Vector& dx,
                                  Vector& r,
                                  const Vector& f,
                                  NLInfo& nlinfo)
{
  const CGInfo& linfo = nlinfo.GetLinearInfo();
  bool lex = linfo.control().status() == "exploded";

  double nn = NewtonNorm(dx);
  double nr = GetSolver(ComputeLevel)->Norm(r);

  if (nn > 1.e30)
    lex = 1;
  if (!(nn >= 0.))
    lex = 1;
  if (nr > 1.e30)
    lex = 1;
  if (!(nr >= 0.))
    lex = 1;

  if (lex) {
    nlinfo.control().status() = "diverged";
    cerr << "linear : " << linfo.control().status() << endl;
    cerr << "nonlinear : " << nn << endl;
    return NewtonNorm(dx);
  }

  double omega = 0.7;
  double relax = 1.;

  GetSolver(ComputeLevel)->SetPeriodicVectorZero(dx);

  GetSolver(ComputeLevel)->Add(x, relax, dx);
  NewtonResidual(r, x, f);
  rr = NewtonNorm(r);

  string message = "";
  for (int iter = 0; iter < nlinfo.user().maxrelax(); iter++) {
    message = nlinfo.check_damping(iter, rr);

    if (message == "ok")
      break;
    if (message == "continue") {
      GetSolver(ComputeLevel)->Add(x, relax * (omega - 1.), dx);

      NewtonResidual(r, x, f);
      rr = NewtonNorm(r);
      relax *= omega;
      continue;
    }
    if (message == "exploded") {
      GetSolver(ComputeLevel)->Add(x, -relax, dx);
      relax = 0.;
      cout << "Damping exploded !!!!!" << endl;
      nlinfo.control().status() = "diverged";
      break;
    }
  }

  // NewtonUpdateShowCompResiduals(nlinfo.control().iteration(), x, r, f,dx);

  return NewtonNorm(dx);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::AssembleMatrix(Matrix& A, Vector& u)
{
  SolutionTransfer(u);
  for (int l = 0; l <= ComputeLevel; l++) {
    GetSolver(l)->MatrixZero(A);
    GetSolver(l)->AssembleMatrix(A, u, 1.);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::AssembleMatrix(Matrix& A, Vector& u, NLInfo& nlinfo)
{
  AssembleMatrix(A, u);
  nlinfo.control().matrixmustbebuild() = 0;
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::ComputeIlu(Matrix& A, Vector& u)
{
  SolutionTransfer(u);
  for (int l = 0; l <= ComputeLevel; l++) {
    GetSolver(l)->ComputeIlu(A, u);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::BoundaryInit(Vector& u) const
{
  for (int l = 0; l < GetSolverPointers().size(); l++)
    GetSolver(l)->BoundaryInit(u);
}

/*-------------------------------------------------------------*/

string
StdMultiLevelSolver::LinearSolve(int level,
                                 const Matrix& A,
                                 Vector& u,
                                 const Vector& b,
                                 CGInfo& info)
{
  ComputeLevel = level;

  GetSolver(ComputeLevel)->HNAverage(u);

  info.reset();

  int clevel = std::max(DataP.CoarseLevel(), 0);
  if (DataP.CoarseLevel() == -1)
    clevel = FinestLevel();

  LinearMg(ComputeLevel, clevel, A, u, b, info);

  GetSolver(ComputeLevel)->SubtractMean(u);

  GetSolver(ComputeLevel)->HNZero(u);
  string status = info.control().status();

  return status;
}

/*-------------------------------------------------------------*/

string
StdMultiLevelSolver::Solve(int level,
                           Matrix& A,
                           Vector& u,
                           const Vector& b,
                           NLInfo& nlinfo)
{
  ComputeLevel = level;
  string status;

  Vector res("_res_");
  Vector cor("_cor_");
  ReInitVector(res);
  ReInitVector(cor);

  assert(DataP.NonLinearSolve() == "newton");
  GetSolver(ComputeLevel)->HNAverage(u);
  newton(A, u, b, res, cor, nlinfo);
  GetSolver(ComputeLevel)->HNZero(u);

  DeleteVector(res);
  DeleteVector(cor);

  return nlinfo.CheckMatrix();
}

/*-------------------------------------------------------------*/

double
StdMultiLevelSolver::ComputeFunctional(Vector& f,
                                       const Vector& u,
                                       const std::string& label)
{
  const Functional* FP = GetFunctionalContainer()->GetFunctional(label);
  return GetSolver(ComputeLevel)->ComputeFunctional(f, u, FP);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::Transfer(int high, int low, Vector& u) const
{
  if (DataP.LinearSolve() == "direct")
    return; // No transfer for direct solver

  for (int l = high; l >= low; l--) {
    GetSolver(l)->HNAverage(u);

    assert(_Interpolator[l - 1]);
    _Interpolator[l - 1]->SolutionTransfer(GetSolver(l - 1)->GetGV(u),
                                           GetSolver(l)->GetGV(u));

    GetSolver(l)->HNZero(u);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::Transfer(Vector& u) const
{
  Transfer(ComputeLevel, 1, u);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::SolutionTransfer(int high, int low, Vector& u) const
{
  if (DataP.LinearSolve() == "direct")
    return; // No transfer for direct solver

  Transfer(high, low, u);

  for (int l = high; l >= low; l--)
    GetSolver(l - 1)->SetBoundaryVector(u);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::SolutionTransfer(Vector& u) const
{
  SolutionTransfer(ComputeLevel, 1, u);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::AssembleDualMatrix(Matrix& A, Vector& u)
{
  for (int l = 0; l < nlevels(); l++) {
    GetSolver(l)->AssembleDualMatrix(A, u, 1.);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::DeleteVector(Vector& v)
{
  for (int l = 0; l < nlevels(); ++l) {
    GetSolver(l)->DeleteVector(v);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::precondition(const Matrix& A, Vector& x, Vector& y)
{
  CGInfo& precinfo = DataP.GetPrecInfo();
  precinfo.reset();
  precinfo.check(0., 0.);

  int clevel = std::max(DataP.CoarseLevel(), 0);
  if (DataP.CoarseLevel() == -1)
    clevel = FinestLevel();

  LinearMg(ComputeLevel, clevel, A, x, y, precinfo);
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::Equ(Vector& dst, double s, const Vector& src) const
{
  for (int l = 0; l < nlevels(); l++) {
    GetSolver(l)->Equ(dst, s, src);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::Zero(Vector& dst) const
{
  for (int l = 0; l < nlevels(); l++) {
    GetSolver(l)->Zero(dst);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::AddNodeVector(const string& name, Vector& gq)
{
  Transfer(ComputeLevel, 1, gq);
  for (int l = 0; l < nlevels(); l++) {
    GetSolver(l)->AddNodeVector(name, gq);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::DeleteNodeVector(const string& name)
{
  for (int l = 0; l < nlevels(); l++) {
    GetSolver(l)->DeleteNodeVector(name);
  }
}

/*-------------------------------------------------------------*/

void
StdMultiLevelSolver::InterpolateCellSolution(Vector& u,
                                             const GlobalVector& uold) const
{
  assert(u.GetType() == "cell");
  GlobalVector& uu = GetSolver()->GetGV(u);

  if (!GetMeshAgent()->Goc2nc()) {
    cerr << "No cell interpolation activated" << endl;
    abort();
  }
  int cells = GetMeshAgent()->ncells();
  uu.ReInit(uold.ncomp(), cells);

  uu.zero();
  // Jetzt Interpolieren wir die Loesung auf das neue Gitter
  for (int i = 0; i < uold.n(); i++) {
    set<int> kinder = GetMeshAgent()->Cello2n(i);
    if (!kinder.empty()) {
      // es wurde nicht vergroebert
      for (set<int>::iterator p = kinder.begin(); p != kinder.end(); p++) {
        for (int c = 0; c < uold.ncomp(); ++c) {
          uu(*p, c) = uold(i, c);
        }
      }
    } else {
      // Es wurde vergroebert
      int patchsize = 1 << GetMeshAgent()->GetMesh(FinestLevel())->dimension();
      uu[GetMeshAgent()->Cello2nFather(i)] +=
        uold[i] / static_cast<double>(patchsize);
    }
  }
}

/*-------------------------------------------------------------*/
} // namespace Gascoigne
