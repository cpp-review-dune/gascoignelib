/**
 *
 * Copyright (C) 2004, 2005, 2006, 2008, 2009, 2011 by the Gascoigne 3D authors
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

#include "basicloop.h"
#include "adaptordata.h"
#include "backup.h"
#include "compose_name.h"
#include "filescanner.h"
#include "gostream.h"
#include "meshagent.h"
#include "monitoring.h"
#include "stdmultilevelsolver.h"
#include <iomanip>

using namespace std;

/*-----------------------------------------*/

double __TIME;
namespace Gascoigne {
extern Timer GlobalTimer;

BasicLoop::BasicLoop()
  : _MA(NULL)
  , _ML(NULL)
  , _SI(NULL)
  , _iter(0)
{
  _reload = "none";
}

BasicLoop::BasicLoop(const ParamFile& paramfile,
                     const ProblemContainer* PC,
                     const FunctionalContainer* FC)
  : _paramfile(paramfile)
{
  //_paramfile = paramfile;
  string s_copy_param_file;
  DataFormatHandler DFH;
  DFH.insert("niter", &_niter, 1);
  DFH.insert("initial", &_initial, "boundary");
  DFH.insert("reload", &_reload, "none");
  DFH.insert("writevtk", &_writeVtk, 1);
  DFH.insert("writebupgup", &_writeBupGup, 1);
  DFH.insert("writeinp", &_writeInp, 0);
  DFH.insert("resultsdir", &_s_resultsdir, "Results");
  DFH.insert("copy_param_file", &s_copy_param_file, "no");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile, "Loop");

  assert((_reload == "none") || (_initial == "file"));
  assert((_reload != "none") || (_initial != "file"));

  GetMeshAgentPointer() = new MeshAgent;
  GetMeshAgent()->BasicInit(_paramfile);

  // GetMultiLevelSolverPointer() = new StdMultiLevelSolver;

  _ML = new StdMultiLevelSolver(GetMeshAgent(), _paramfile, PC, FC);

  GetMultiLevelSolver()->GetMonitor().set_directory(_s_resultsdir);
  GetSolverInfosPointer() = new SolverInfos;
  GetSolverInfos()->BasicInit(_paramfile);
}

/*-----------------------------------------*/

BasicLoop::~BasicLoop()
{
  if (_MA != NULL) {
    delete _MA;
    _MA = NULL;
  }
  if (_ML != NULL) {
    delete _ML;
    _ML = NULL;
  }
  if (_SI != NULL) {
    delete _SI;
    _SI = NULL;
  }
}

/*-----------------------------------------*/

void
BasicLoop::ClockOutput() const
{
  cout << "********************************************************************"
          "****\n\n";
  cout << "BasicLoop\t\tTIME\n";
  cout << "  NewMesh\t\t" << _clock_newmesh.read() << endl;
  cout << "  Solve\t\t\t" << _clock_solve.read() << endl;
  cout << "  Write\t\t\t" << _clock_write.read() << endl;
}

/*-----------------------------------------*/

void
BasicLoop::BasicInit(const ParamFile& paramfile,
                     const ProblemContainer* PC,
                     const FunctionalContainer* FC)
{
  _paramfile = paramfile;
  string s_copy_param_file;

  DataFormatHandler DFH;
  DFH.insert("niter", &_niter, 1);
  DFH.insert("initial", &_initial, "boundary");
  DFH.insert("reload", &_reload, "none");
  DFH.insert("writevtk", &_writeVtk, 1);
  DFH.insert("writebupgup", &_writeBupGup, 1);
  DFH.insert("writeinp", &_writeInp, 0);
  DFH.insert("resultsdir", &_s_resultsdir, "Results");
  DFH.insert("copy_param_file", &s_copy_param_file, "no");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile, "Loop");

  assert((_reload == "none") || (_initial == "file"));

  if ((_reload != "none") && (_initial != "file")) {
    cerr << "Please, add 'initial file' to Block Loop" << endl;
    _initial = "file";
  }
  assert((_reload != "none") || (_initial != "file"));

  if (GetMeshAgentPointer() == NULL) {
    GetMeshAgentPointer() = new MeshAgent;
  }
  GetMeshAgent()->BasicInit(_paramfile);

  if (GetMultiLevelSolverPointer() == NULL) {
    GetMultiLevelSolverPointer() = new StdMultiLevelSolver;
  }
  assert(GetMultiLevelSolver());

  GetMultiLevelSolver()->BasicInit(GetMeshAgent(), _paramfile, PC, FC);

  if (GetSolverInfosPointer() == NULL) {
    GetSolverInfosPointer() = new SolverInfos;
  }
  assert(GetSolverInfos());
  GetSolverInfos()->BasicInit(_paramfile);
}

/*-------------------------------------------------------*/

void
BasicLoop::PrintMeshInformation(int outputlevel) const
{
  cout << " [l,n,c] " << GetMeshAgent()->nlevels() << " "
       << GetMeshAgent()->nnodes();
  cout << " " << GetMeshAgent()->ncells() << endl;

  if (outputlevel) {
    for (int l = 0; l < GetMeshAgent()->nlevels(); l++) {
      const GascoigneMesh* M = GetMeshAgent()->GetMesh(l);
      cout << l << " [n,c] " << M->nnodes() << " " << M->ncells() << endl;
    }
  }
}

/*-------------------------------------------------------*/

void
BasicLoop::Output(const Vector& u, string name) const
{
  if (_writeVtk > 0) {
    if (_iter % _writeVtk == 0)
      GetMultiLevelSolver()->GetSolver()->Visu(name, u, _iter);
  }
  if (_writeBupGup > 0) {
    if (_iter % _writeBupGup == 0)
      WriteMeshAndSolution(name, u);
  }
  if (_writeInp > 0) {
    if (_iter % _writeInp == 0)
      WriteMeshInp(name);
  }
}

/*-------------------------------------------------*/

void
BasicLoop::WriteMeshAndSolution(const string& filename, const Vector& u) const
{
  string name;
  name = filename;
  //   name = filename + "_value";
  compose_name(name, _iter);
  GetMultiLevelSolver()->GetSolver()->Write(u, name);
  cout << "[" << name << ".bup]";

  //   name = filename + "_mesh";
  //   compose_name(name,_iter);
  GetMeshAgent()->write_gup(name);
  cout << " [" << name << ".gup]" << endl;
}

/*-------------------------------------------------*/

void
BasicLoop::WriteSolution(const Vector& u) const
{
  _clock_write.start();
  string filename = _s_resultsdir + "/solution";
  compose_name(filename, _iter);
  GetMultiLevelSolver()->GetSolver()->Write(u, filename);
  cout << "[" << filename << ".bup]" << endl;
  _clock_write.stop();
}

/*-------------------------------------------------*/

void
BasicLoop::WriteMesh() const
{
  _clock_write.start();
  string filename = _s_resultsdir + "/mesh";
  compose_name(filename, _iter);
  GetMeshAgent()->write_gup(filename);
  cout << " [" << filename << ".gup]" << endl;
  _clock_write.stop();
}

/*-------------------------------------------------*/

void
BasicLoop::WriteMeshInp(const string& name) const
{
  _clock_write.start();
  string filename = name;
  compose_name(filename, _iter);
  GetMeshAgent()->write_inp(filename);
  cout << " [" << filename << ".inp]" << endl;
  _clock_write.stop();
}

/*-------------------------------------------------*/

void
BasicLoop::InitSolution(Vector& u)
{
  GetMultiLevelSolver()->GetSolver()->Zero(u);

  if (_initial == "analytic")
    GetMultiLevelSolver()->GetSolver()->SolutionInit(u);
  else if (_initial == "file")
    GetMultiLevelSolver()->GetSolver()->Read(u, _reload);
  else if (_initial == "boundary")
    GetMultiLevelSolver()->GetSolver()->BoundaryInit(u);
  else if (_initial != "zero") {
    cerr << "BasicLoop::InitSolution():\npossible values for init: "
            "\"analytic\", \"file\", "
            "\"boundary\", \"zero\"\n";
    abort();
  }
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  GetMultiLevelSolver()->GetSolver()->SubtractMean(u);
  GetMultiLevelSolver()->GetSolver()->Visu(_s_resultsdir + "/solve", u, 0);
}

/*-------------------------------------------------*/

string
BasicLoop::Solve(Matrix& A, Vector& u, Vector& f, string name)
{
  _clock_solve.start();

  GetMultiLevelSolver()->GetSolver()->Zero(f);
  GetMultiLevelSolver()->GetSolver()->Rhs(f);

  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  // set offset first, so nodes that are both periodic and dirichlet will become
  // dirichlet

  GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);

  string status =
    GetMultiLevelSolver()->Solve(A, u, f, GetSolverInfos()->GetNLInfo());
  _clock_solve.stop();

  _clock_write.start();
  Output(u, name);
  _clock_write.stop();

  return status;
}

/*-------------------------------------------------*/

void
BasicLoop::ComputeGlobalErrors(const Vector& u)
{
  GetMultiLevelSolver()->GetSolver()->ComputeError(u, _GlobalErr);
  if (_GlobalErr.size() > 0) {
    cout.precision(6);
    for (int c = 0; c < _GlobalErr.ncomp(); c++) {
      cout << "\nGlobalErrors [" << c << "]  l2,h1,l8 ";
      for (int i = 0; i < _GlobalErr.n(); i++)
        cout << " " << _GlobalErr(i, c);
    }
    cout << endl;
  }
}

/*-------------------------------------------------------*/

void
BasicLoop::CopyVector(GlobalVector& dst, Vector& src)
{
  GetMultiLevelSolver()->GetSolver()->HNAverage(src);

  int nn = GetMultiLevelSolver()->GetSolver()->GetGV(src).n();
  int cc = GetMultiLevelSolver()->GetSolver()->GetGV(src).ncomp();

  dst.ncomp() = cc;
  dst.resize(nn);

  dst.equ(1., GetMultiLevelSolver()->GetSolver()->GetGV(src));

  GetMultiLevelSolver()->GetSolver()->HNZero(src);
}

/*-------------------------------------------------*/

void
BasicLoop::CopyVector(Vector& dst, GlobalVector& src)
{
  int nn = src.n();
  int cc = src.ncomp();

  GetMultiLevelSolver()->GetSolver()->GetGV(dst).ncomp() = cc;
  GetMultiLevelSolver()->GetSolver()->GetGV(dst).resize(nn);
  GetMultiLevelSolver()->GetSolver()->GetGV(dst).equ(1., src);
}

/*-------------------------------------------------*/

void
BasicLoop::run(const std::string& problemlabel)
{
  Matrix A("A");
  Vector u("u"), f("f");
  GlobalVector ualt;

  Monitoring Moning;

  for (_iter = 1; _iter <= _niter; _iter++) {
    cout << "\n================== " << _iter << " ================";
    PrintMeshInformation();
    Moning.SetMeshInformation(
      _iter, GetMeshAgent()->nnodes(), GetMeshAgent()->ncells());

    _clock_newmesh.start();

    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    GetMultiLevelSolver()->ReInit();
    GetMultiLevelSolver()->SetProblem(problemlabel);
    GetMultiLevelSolver()->ReInitMatrix(A);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(f);
    GetMultiLevelSolver()->InterpolateSolution(u, ualt);
    GetMultiLevelSolver()->GetSolver()->Visu(
      _s_resultsdir + "/interpolate", u, _iter);

    _clock_newmesh.stop();

    if (_iter == 1) {
      GetMultiLevelSolver()->GetSolver()->OutputSettings();
      InitSolution(u);
    }

    Solve(A, u, f);
    ComputeGlobalErrors(u);

    if (_iter < _niter) {
      CopyVector(ualt, u);

      GetMeshAgent()->global_refine(1);
    }
  }
}

void
BasicLoop::timerun(const std::string& problemlabel)
{
  Vector u("u");
  Vector f("f");
  Vector old("old");
  Matrix A("A");

  GetMultiLevelSolver()->ReInit();
  GetMultiLevelSolver()->SetProblem(problemlabel);
  GetMultiLevelSolver()->ReInitMatrix(A);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(f);
  GetMultiLevelSolver()->ReInitVector(old);
  InitSolution(u);

  double dt;
  DataFormatHandler DFH;
  DFH.insert("dt", &dt, 1.);
  FileScanner FS(DFH, _paramfile, "Equation");
  assert(dt > 0);

  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

  PrintMeshInformation();
  GlobalTimer.reset();

  __TIME = 0;
  for (_iter = 1; _iter <= _niter; _iter++) {
    GlobalTimer.start("iteration");

    std::cout << "\n******************** Time Step " << _iter << "\t" << __TIME
              << " -> " << __TIME + dt << std::endl;
    __TIME += dt;
    GetMultiLevelSolver()->Equ(old, 1., u);

    GetMultiLevelSolver()->AddNodeVector("OLD", old);
    Solve(A, u, f);
    GetMultiLevelSolver()->DeleteNodeVector("OLD");
    GlobalTimer.stop("iteration");
    // GlobalTimer.print();
  }
  GlobalTimer.print();
}
} // namespace Gascoigne

/*-------------------------------------------------*/
