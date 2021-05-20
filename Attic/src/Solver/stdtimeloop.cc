/**
 *
 * Copyright (C) 2004, 2005, 2006, 2009 by the Gascoigne 3D authors
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

#include "stdtimeloop.h"
#include "filescanner.h"
#include "stdtimesolver.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {
void
StdTimeLoop::BasicInit(const ParamFile* paramfile,
                       const ProblemContainer* PC,
                       const FunctionalContainer* FC)
{
  StdLoop::BasicInit(paramfile, PC, FC);

  double tbegin, tend, deltat, theta;
  int neuler;
  string scheme;

  DataFormatHandler DFH;
  DFH.insert("dt", &deltat, 1.);
  DFH.insert("tbegin", &tbegin, 0.);
  DFH.insert("tend", &tend, 1.e4);
  DFH.insert("neuler", &neuler, 10);
  DFH.insert("scheme", &scheme, "Euler");
  DFH.insert("theta", &theta, 0.5);
  DFH.insert("reload", &_reload);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile, "Loop");

  _timeinfo.ReInit(tbegin, tend, deltat, scheme, neuler, theta);
}

/*-------------------------------------------------*/

string
StdTimeLoop::SolveTimePrimal(Matrix& A, Vector& u, Vector& f)
{
  GetTimeSolver()->SetBoundaryVector(f);
  GetTimeSolver()->SetPeriodicVector(u);
  GetTimeSolver()->SetBoundaryVector(u);
  string status =
    GetMultiLevelSolver()->Solve(A, u, f, GetSolverInfos()->GetNLInfo());

  return status;
}

/*-------------------------------------------------*/

void
StdTimeLoop::adaptive_run(const std::string& problemlabel)
{
  Matrix A("A");
  Vector u("u"), f("f");
  GlobalVector ualt;

  DoubleVector eta;

  for (_iter = 1; _iter <= _niter; _iter++) {
    GetMultiLevelSolver()->ReInit(problemlabel);
    GetMultiLevelSolver()->ReInitMatrix(A);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(f);
    TimeInfoBroadcast();

    GetMultiLevelSolver()->InterpolateSolution(u, ualt);

    if (_iter == 1) {
      GetTimeSolver()->OutputSettings();
      InitSolution(u);
    }

    cout << "\n================== " << _iter << "================";
    cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " "
         << GetMeshAgent()->nnodes();
    cout << " " << GetMeshAgent()->ncells() << endl;

    // umschalten von Euler ?
    //
    _timeinfo.SpecifyScheme(_iter);
    TimeInfoBroadcast();
    //
    // rhs fuer alten Zeitschritt
    //
    GetTimeSolver()->GetGV(f).zero();
    GetTimeSolver()->TimeRhsOperator(f, u);
    GetTimeSolver()->TimeRhs(1, f);

    // neuer Zeitschritt
    //
    _timeinfo.iteration(_iter);
    TimeInfoBroadcast();

    GetTimeSolver()->TimeRhs(2, f);
    SolveTimePrimal(A, u, f);
    Output(u, _s_resultsdir + "/solve");

    StdSolver* S = dynamic_cast<StdSolver*>(GetTimeSolver());
    assert(S);
    if (_iter < _niter) {
      CopyVector(ualt, u);

      AdaptMesh(eta);
    }
  }
}

/*-------------------------------------------------*/

void
StdTimeLoop::TimeInfoBroadcast()
{
  for (int l = 0; l < GetMultiLevelSolver()->nlevels(); l++) {
    StdTimeSolver* TS =
      dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver(l));
    assert(TS);
    TS->SetTimeData(_timeinfo.dt(),
                    _timeinfo.theta(),
                    _timeinfo.time(),
                    _timeinfo.oldrhs(),
                    _timeinfo.rhs());
  }
}

/*-------------------------------------------------*/

void
StdTimeLoop::InitSolution(Vector& u)
{
  if (_initial == "analytic") {
    Vector f("ff");
    GetMultiLevelSolver()->ReInitVector(f);
    GetTimeSolver()->L2Projection(u, f);
    GetMultiLevelSolver()->DeleteVector(f);
    GetTimeSolver()->Write(u, _s_resultsdir + "/initialu");
  } else {
    StdLoop::InitSolution(u);
  }
}

/*-------------------------------------------------*/

void
StdTimeLoop::run(const std::string& problemlabel)
{
  Vector u("u"), f("f");
  Matrix A("A");

  DoubleVector eta;

  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitMatrix(A);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(f);

  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " "
       << GetMeshAgent()->ncells() << endl;
  GetTimeSolver()->OutputSettings();

  TimeInfoBroadcast();

  // Anfangswerte
  InitSolution(u);

  GetTimeSolver()->SetPeriodicVector(u);
  GetTimeSolver()->SetBoundaryVector(u);
  GetTimeSolver()->Visu(_s_resultsdir + "/solve", u, 0);

  for (_iter = 1; _iter <= _niter; _iter++) {
    // umschalten von Euler ?
    //
    _timeinfo.SpecifyScheme(_iter);
    TimeInfoBroadcast();
    //
    // rhs fuer alten Zeitschritt
    //
    GetTimeSolver()->GetGV(f).zero();
    GetTimeSolver()->TimeRhsOperator(f, u);
    GetTimeSolver()->TimeRhs(1, f);

    // neuer Zeitschritt
    //
    _timeinfo.iteration(_iter);
    TimeInfoBroadcast();

    GetTimeSolver()->TimeRhs(2, f);

    SolveTimePrimal(A, u, f);
    Output(u, _s_resultsdir + "/solve");

    Functionals(u, f);
  }
}
} // namespace Gascoigne
