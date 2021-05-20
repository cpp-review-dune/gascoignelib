/**
 *
 * Copyright (C) 2008, 2010 by the Gascoigne 3D authors
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

#include "onelevelalgorithm.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {

/*-----------------------------------------*/

void
OneLevelAlgorithm::BasicInit(const ParamFile* paramfile,
                             const NumericInterface* NI,
                             const ProblemContainer* PC)
{
  Algorithm::BasicInit(paramfile, NI);

  _S = GetNumeric()->NewSolver();

  GetSolver()->SetDiscretization(*GetNumeric()->NewDiscretization());
  GetSolver()->BasicInit(paramfile, GetMeshAgent()->GetDimension());

  _PC = PC;
}

/*-----------------------------------------*/

void
OneLevelAlgorithm::Precondition(Vector& x, Vector& y)
{
  CGInfo pinfo;
  pinfo.user().tol() = 1.e-12;
  pinfo.user().globaltol() = 1.e-12;
  pinfo.user().maxiter() = 1;
  pinfo.user().printstep() = 0;
  pinfo.user().text() = "PrecInfo";

  JacobiSolver(x, y, pinfo);
}

/*-----------------------------------------*/

void
OneLevelAlgorithm::IluSolver(Vector& du, const Vector& f, CGInfo& info)
{
  Vector help("help");
  ReInitVector(help);

  bool ok = info.check(GetSolver()->Norm(f), 0.);
  for (int iter = 0; !ok; iter++) {
    GetSolver()->MatrixResidual(help, du, f);
    double rnorm = GetSolver()->Norm(help);

    StdSolver* SS = dynamic_cast<StdSolver*>(GetSolver());
    SS->GetIlu()->solve(GetSolver()->GetGV(help));
    double cnorm = GetSolver()->Norm(help);

    GetSolver()->Add(du, 1., help);
    ok = info.check(rnorm, cnorm);
  }
  DeleteVector(help);
}

/*-----------------------------------------*/

void
OneLevelAlgorithm::JacobiSolver(Vector& du, const Vector& f, CGInfo& info)
{
  Vector help("help");
  ReInitVector(help);

  bool ok = info.check(GetSolver()->Norm(f), 0.);
  for (int iter = 0; !ok; iter++) {
    GetSolver()->MatrixResidual(help, du, f);
    double rnorm = GetSolver()->Norm(help);

    StdSolver* SS = dynamic_cast<StdSolver*>(GetSolver());
    SS->GetMatrix()->Jacobi(GetSolver()->GetGV(help));

    double cnorm = GetSolver()->Norm(help);

    GetSolver()->Add(du, 1., help);
    ok = info.check(rnorm, cnorm);
  }
  DeleteVector(help);
}

/*-----------------------------------------*/

void
OneLevelAlgorithm::RunLinear(const std::string& problemlabel)
{
  GetSolver()->NewMesh(GetMeshAgent()->GetMesh(0));

  GetSolver()->SetProblem(*_PC->GetProblem(problemlabel));

  GetSolver()->OutputSettings();
  PrintMeshInformation();

  Vector u("u"), f("f"), du("du");

  ReInitVector(u);
  ReInitVector(f);
  ReInitVector(du);

  GetSolver()->Zero(u);
  GetSolver()->SetBoundaryVector(u);

  // Compute RHS

  GetSolver()->Zero(f);
  GetSolver()->Rhs(f);
  GetSolver()->SetBoundaryVector(f);

  // Compute Residual

  GetSolver()->HNAverage(u);
  GetSolver()->Form(f, u, -1.);
  GetSolver()->SetBoundaryVectorZero(f);

  // Assemble Matrix and ILU

  GetSolver()->RegisterMatrix();
  GetSolver()->ReInitMatrix();
  GetSolver()->MatrixZero();
  GetSolver()->AssembleMatrix(u, 1.);
  GetSolver()->ComputeIlu();

  // Solve Linear System

  GetSolver()->Zero(du);

  CGInfo& info = GetSolverInfos()->GetLInfo();

  IluSolver(du, f, info);

  cout << endl << "Linear solver " << info.control().status() << endl << endl;
  GetSolver()->Add(u, 1., du);
  GetSolver()->HNZero(u);

  GetSolver()->Visu("Results/onelevel", u, 0);

  DeleteVector(u);
  DeleteVector(f);
  DeleteVector(du);
}

/*-----------------------------------------*/

void
OneLevelAlgorithm::AssembleMatrixAndIlu(Vector& u)
{
  GetSolver()->MatrixZero();
  GetSolver()->AssembleMatrix(u, 1.);
  GetSolver()->ComputeIlu(u);
}

/*-----------------------------------------*/

void
OneLevelAlgorithm::LinearSolve(Vector& du, const Vector& y, CGInfo& cginfo)
{
  cginfo.reset();
  IluSolver(du, y, cginfo);
}

/*-----------------------------------------*/

void
OneLevelAlgorithm::RunNonLinear(const std::string& problemlabel)
{
  GetSolver()->NewMesh(GetMeshAgent()->GetMesh(0));

  GetSolver()->SetProblem(*_PC->GetProblem(problemlabel));

  GetSolver()->OutputSettings();
  PrintMeshInformation();

  Vector u("u"), f("f");

  ReInitVector(u);
  ReInitVector(f);

  GetSolver()->Zero(u);
  GetSolver()->SolutionInit(u);
  GetSolver()->SubtractMean(u);
  GetSolver()->SetBoundaryVector(u);

  // Compute RHS

  GetSolver()->Zero(f);
  GetSolver()->Rhs(f);
  GetSolver()->SetBoundaryVector(f);

  // Memory Matrix

  GetSolver()->RegisterMatrix();
  GetSolver()->ReInitMatrix();
  GetSolver()->MatrixZero();

  NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();

  Newton(u, f, nlinfo);

  GetSolver()->Visu("Results/onelevel", u, 0);

  DeleteVector(u);
  DeleteVector(f);
}
} // namespace Gascoigne
