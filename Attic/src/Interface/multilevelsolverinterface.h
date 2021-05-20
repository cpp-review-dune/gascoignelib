/**
 *
 * Copyright (C) 2004, 2005, 2006, 2008 by the Gascoigne 3D authors
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

#ifndef __MultiLevelSolverInterface_h
#define __MultiLevelSolverInterface_h

#include "functionalcontainer.h"
#include "meshagentinterface.h"
#include "monitor.h"
#include "nlinfo.h"
#include "paramfile.h"
#include "problemcontainer.h"
#include "solverinterface.h"
#include "vectorinterface.h"

namespace Gascoigne {

/////////////////////////////////////////////
///
///@brief
///  ... comments MultiLevelSolverInterface

///
///
/////////////////////////////////////////////

class MultiLevelSolverInterface
{
private:
protected:
public:
  MultiLevelSolverInterface() {}
  virtual ~MultiLevelSolverInterface() {}

  virtual std::string GetName() const = 0;
  virtual void BasicInit(const MeshAgentInterface* GMGM,
                         const ParamFile* paramfile,
                         const ProblemContainer* PC,
                         const FunctionalContainer* FC = NULL) = 0;
  virtual void SetProblem(const std::string& problemlabel) = 0;
  virtual void ReInit(const std::string& problemlabel) = 0;
  virtual void SetMonitorPtr(Monitor* mon) = 0;

  virtual const DoubleVector GetExactValues() const = 0;
  virtual const DoubleVector ComputeFunctionals(Vector& f, const Vector& u) = 0;
  virtual const DoubleVector ComputeFunctionals(Vector& f,
                                                const Vector& u,
                                                FunctionalContainer* FC) = 0;

  virtual void ReInitMatrix() = 0;

  virtual int nlevels() const = 0;

  virtual SolverInterface* GetSolver(int l) = 0;
  virtual const SolverInterface* GetSolver(int l) const = 0;
  virtual SolverInterface* GetSolver() = 0;
  virtual const SolverInterface* GetSolver() const = 0;
  virtual const ProblemContainer* GetProblemContainer() const = 0;

  //      virtual void SetState(const std::string& s)=0;
  virtual void AssembleMatrix(Vector& u, NLInfo& nlinfo) = 0;
  virtual void AssembleMatrix(Vector& u) = 0;
  virtual void ComputeIlu(Vector& u) = 0;
  virtual void ComputeIlu() = 0;

  virtual void BoundaryInit(Vector& u) const = 0;

  //
  /// vector - manamgement
  //

  virtual void DeleteVector(Vector& g) = 0;
  virtual void RegisterVectors() = 0;
  virtual void RegisterMatrix() = 0;
  virtual void ReInitVector(Vector& v) = 0;
  virtual void ReInitVector(Vector& v, int comp) = 0;

  //
  /// vector
  //

  virtual std::string LinearSolve(int level,
                                  Vector& u,
                                  const Vector& b,
                                  CGInfo& info) = 0;
  virtual std::string LinearSolve(Vector& u, const Vector& b, CGInfo& info)
  {
    return LinearSolve(nlevels() - 1, u, b, info);
  }

  virtual std::string Solve(int level,
                            Vector& x,
                            const Vector& b,
                            NLInfo& nlinfo) = 0;
  virtual std::string Solve(Vector& x, const Vector& b, NLInfo& nlinfo)
  {
    return Solve(nlevels() - 1, x, b, nlinfo);
  }
  virtual void InterpolateSolution(Vector& u,
                                   const GlobalVector& uold) const = 0;
  virtual void InterpolateCellSolution(Vector& u,
                                       const GlobalVector& uold) const = 0;
  virtual double ComputeFunctional(Vector& f,
                                   const Vector& u,
                                   const std::string& label) = 0;

  virtual void AssembleDualMatrix(Vector& u) = 0;
  virtual void vmulteq(Vector& y, const Vector& x) const = 0;

  virtual void Equ(Vector& dst, double s, const Vector& src) const = 0;
  virtual void Zero(Vector& dst) const = 0;
  virtual void AddNodeVector(const std::string& name, Vector& q) = 0;
  virtual void DeleteNodeVector(const std::string& q) = 0;

  virtual void SolutionTransfer(Vector& u) const = 0;
  virtual void Transfer(Vector& u) const = 0;

  virtual void newton(Vector& u,
                      const Vector& f,
                      Vector& r,
                      Vector& w,
                      NLInfo& info) = 0;

  virtual std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers() = 0;
  virtual const std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers()
    const = 0;
};
} // namespace Gascoigne

#endif
