/**
 *
 * Copyright (C) 2004, 2005, 2007 by the Gascoigne 3D authors
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

#ifndef __BasicLoop_h
#define __BasicLoop_h

#include "../Common/extrapolator.h"
#include "../Common/paramfile.h"
#include "../Common/stopwatch.h"
#include "../Discretization/Q1/visualization.h"
#include "../Discretization/Q1/visudatacompvector.h"
#include "../Discretization/Q1/visudatanvector.h"
#include "../Interface/matrix.h"
#include "../Interface/problemdescriptorinterface.h"
#include "../Interface/vectorinterface.h"
#include "../Mesh/adaptordata.h"
#include "../Mesh/meshagent.h"
#include "../Problems/functionalcontainer.h"

#include "monitor.h"
#include "solverinfos.h"
#include "stdmultilevelsolver.h"

/*-----------------------------------------*/

namespace Gascoigne {
//////////////////////////////////////////////
//
///@brief
/// Implementation of diverse outer loops: mesh adaptation, time stepping...

///
///
//////////////////////////////////////////////

class BasicLoop
{
private:
  MeshAgent* _MA;
  StdMultiLevelSolver* _ML;
  SolverInfos* _SI;

protected:
  void WriteMeshAndSolution(const std::string& filename, const Vector& u) const;
  void WriteSolution(const Vector& u) const;
  void WriteMesh() const;
  void WriteMeshInp(const std::string& name) const;

  virtual MeshAgent*& GetMeshAgentPointer() { return _MA; }
  virtual StdMultiLevelSolver*& GetMultiLevelSolverPointer() { return _ML; }

  virtual const MeshAgent* GetMeshAgent() const
  {
    assert(_MA);
    return _MA;
  }
  virtual const StdMultiLevelSolver* GetMultiLevelSolver() const
  {
    assert(_ML);
    return _ML;
  }

  virtual MeshAgent* GetMeshAgent()
  {
    assert(_MA);
    return _MA;
  }
  virtual StdMultiLevelSolver* GetMultiLevelSolver()
  {
    assert(_ML);
    return _ML;
  }

  virtual SolverInfos*& GetSolverInfosPointer() { return _SI; }
  virtual SolverInfos* GetSolverInfos()
  {
    assert(_SI);
    return _SI;
  }
  virtual const SolverInfos* GetSolverInfos() const
  {
    assert(_SI);
    return _SI;
  }

  mutable StopWatch _clock_newmesh, _clock_solve, _clock_functionals,
    _clock_write, _clock_estimate;

  int _niter, _iter;
  int _writeVtk;
  int _writeBupGup;
  int _writeInp;

  std::string _reload, _initial;
  std::string _s_resultsdir;
  ParamFile _paramfile;

  GlobalVector _GlobalErr;

  // new vectors

  virtual std::string Solve(Matrix& A, Vector& u, Vector& f, std::string name);
  virtual std::string Solve(Matrix& A, Vector& u, Vector& f)
  {
    return Solve(A, u, f, _s_resultsdir + "/solve");
  }

  virtual void PrintMeshInformation(int outputlevel = 0) const;

  virtual void Output(const Vector& u, std::string name) const;
  virtual void Output(const Vector& u) { Output(u, _s_resultsdir + "/solve"); }

  virtual void ComputeGlobalErrors(const Vector& u);

  virtual void InitSolution(Vector& u);
  virtual void CopyVector(GlobalVector& dst, Vector& src);
  virtual void CopyVector(Vector& dst, GlobalVector& src);

public:
  BasicLoop();
  BasicLoop(const ParamFile& paramfile,
            const ProblemContainer* PC,
            const FunctionalContainer* FC);

  virtual ~BasicLoop();

  virtual void BasicInit(const ParamFile& paramfile,
                         const ProblemContainer* PC,
                         const FunctionalContainer* FC);

  void run(const std::string& problemlabel);
  void timerun(const std::string& problemlabel);
  void ClockOutput() const;
};
} // namespace Gascoigne

/*-----------------------------------------*/

#endif
