/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 by the Gascoigne 3D authors
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

#ifndef __StdMultiLevelSolver_h
#define __StdMultiLevelSolver_h

#include "functionalcontainer.h"
#include "meshagent.h"
#include "mginterpolatorinterface.h"
#include "monitor.h"
#include "problemcontainer.h"
#include "problemdescriptorinterface.h"
#include "stdmultilevelsolverdata.h"
#include "stdsolver.h"
#include "stopwatch.h"
#include "vectorinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne {
//////////////////////////////////////////////
///
///@brief
/// Default nonlinear MultilevelSolver

/// - stores MultiGridMeshInterace
/// - stores array of MGInterpolator
/// - stores array of SolverInterface
///
//////////////////////////////////////////////

class StdMultiLevelSolver
{
protected: ///!!!!
  std::vector<StdSolver*> _SP;
  const MeshAgent* _MAP;
  std::vector<MgInterpolatorInterface*> _Interpolator;

  mutable StopWatch _clock_residual, _clock_solve;
  mutable IndexType ComputeLevel;
  mutable IndexType oldnlevels;

  ParamFile _paramfile;

  Monitor MON;

  const ProblemDescriptorInterface* _PD;
  const ProblemContainer* _PC;
  const FunctionalContainer* _FC;

  StdMultiLevelSolverData DataP;

protected:
public:
  //////////////////////////////////////////////////
  // Constructor, Init

  StdMultiLevelSolver();
  explicit StdMultiLevelSolver(const MeshAgent* GMGM,
                               const ParamFile& paramfile,
                               const ProblemContainer* PC,
                               const FunctionalContainer* FC = NULL);
  virtual ~StdMultiLevelSolver();

  virtual std::string GetName() const { return "StdMultiLevelSolver"; }

  virtual void BasicInit(const MeshAgent* GMGM,
                         const ParamFile& paramfile,
                         const ProblemContainer* PC,
                         const FunctionalContainer* FC = NULL);

  virtual void ReInitMatrix(const Matrix& A);
  virtual void ReInitVector(Vector& v);
  virtual void ReInitVector(Vector& v, IndexType comp);

  //////////////////////////////////////////////////
  // Access

  virtual const int GetComputeLevel() const { return ComputeLevel; }

  virtual const MeshAgent* GetMeshAgent() const { return _MAP; }
  virtual std::vector<StdSolver*>& GetSolverPointers() { return _SP; }
  virtual const std::vector<StdSolver*>& GetSolverPointers() const
  {
    return _SP;
  }
  virtual const ProblemDescriptorInterface* GetProblemDescriptor() const
  {
    return _PD;
  }

  virtual StdSolver*& GetSolverPointer(IndexType l)
  {
    assert(l < _SP.size());
    return _SP[l];
  }

  //////////////////////////////////////////////////
  // Solver
  virtual void NewSolvers();

  virtual StdSolver* NewSolver(IndexType solverlevel);
  virtual void NewMgInterpolator();
  virtual void SolverNewMesh();
  virtual void SetComputeLevel(IndexType level) { ComputeLevel = level; }
  virtual std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers()
  {
    return _Interpolator;
  }
  virtual const std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers()
    const
  {
    return _Interpolator;
  }
  virtual const ProblemContainer* GetProblemContainer() const
  {
    assert(_PC);
    return _PC;
  }
  virtual void SetProblemContainer(const ProblemContainer* PC) { _PC = PC; }
  virtual const FunctionalContainer* GetFunctionalContainer() const
  {
    return _FC;
  }
  virtual void SetFunctionalContainer(const FunctionalContainer* FC)
  {
    _FC = FC;
  }

  virtual StdSolver* GetSolver(IndexType l)
  {
    assert(l < _SP.size());
    return _SP[l];
  }
  virtual const StdSolver* GetSolver(IndexType l) const
  {
    assert(l < _SP.size());
    return _SP[l];
  }
  virtual StdSolver* GetSolver()
  {
    assert(_SP.size() == nlevels());
    return _SP[FinestLevel()];
  }
  virtual const StdSolver* GetSolver() const
  {
    assert(_SP.size() == nlevels());
    return _SP[FinestLevel()];
  }

  //////////////////////////////////////////////////
  // Functionals
  virtual const std::vector<std::string> GetFunctionalNames() const;
  virtual const DoubleVector GetExactValues() const;
  virtual const DoubleVector ComputeFunctionals(Vector& f, const Vector& u);
  virtual const DoubleVector ComputeFunctionals(Vector& f,
                                                const Vector& u,
                                                FunctionalContainer* FC);

  //////////////////////////////////////////////////
  // Newton and Multigrid
  virtual double NewtonNorm(const Vector& u) const
  {
    return GetSolver(ComputeLevel)->NewtonNorm(u);
  }
  virtual void mgstep(std::vector<double>& res,
                      std::vector<double>& rw,
                      IndexType l,
                      IndexType maxl,
                      IndexType minl,
                      std::string& p0,
                      std::string p,
                      const Matrix& A,
                      Vector& u,
                      Vector& b,
                      Vector& v);

  virtual void Cg(Vector& x, const Vector& f, CGInfo& info);
  virtual void Gmres(const Matrix& A, Vector& x, const Vector& f, CGInfo& info);

  virtual void SolutionTransfer(IndexType high, IndexType low, Vector& u) const;
  virtual void Transfer(IndexType high, IndexType low, Vector& u) const;

  virtual void ViewProtocoll() const;

  virtual IndexType nlevels() const
  {
    assert(GetMeshAgent());
    return GetMeshAgent()->nlevels();
  }
  virtual IndexType FinestLevel() const { return nlevels() - 1; }
  virtual IndexType CoarsestLevel() const { return 0; }

  Monitor& GetMonitor() { return MON; }

  const Monitor& GetMonitor() const { return MON; }

  // ReInit is to be called, whenever the mesh (hierarchy) has changed
  virtual void ReInit();
  // SetProblem announces the current ProblemDescriptor to all solvers in the
  // hierarchy
  virtual void SetProblem(const std::string& label);

  // neue vektoren

  virtual std::string LinearSolve(IndexType level,
                                  const Matrix& A,
                                  Vector& u,
                                  const Vector& b,
                                  CGInfo& info);
  virtual std::string Solve(IndexType level,
                            Matrix& A,
                            Vector& x,
                            const Vector& b,
                            NLInfo& nlinfo);
  virtual std::string Solve(Matrix& A,
                            Vector& x,
                            const Vector& b,
                            NLInfo& nlinfo)
  {
    return Solve(nlevels() - 1, A, x, b, nlinfo);
  }

  virtual void InterpolateSolution(Vector& u, const GlobalVector& uold) const;
  virtual void InterpolateCellSolution(Vector& u,
                                       const GlobalVector& uold) const;

  virtual void NewtonVectorZero(Vector& w) const;
  virtual double NewtonResidual(Vector& y,
                                const Vector& x,
                                const Vector& b) const;
  virtual double NewtonUpdate(double& rr,
                              Vector& x,
                              Vector& dx,
                              Vector& r,
                              const Vector& f,
                              NLInfo& nlinfo);
  virtual void NewtonLinearSolve(const Matrix& A,
                                 Vector& x,
                                 const Vector& b,
                                 CGInfo& info);
  virtual void NewtonMatrixControl(Matrix& A, Vector& u, NLInfo& nlinfo);
  virtual void NewtonOutput(NLInfo& nlinfo) const;
  virtual void NewtonPreProcess(Vector& u, const Vector& f, NLInfo& info) const;
  virtual void NewtonPostProcess(Vector& u,
                                 const Vector& f,
                                 NLInfo& info) const;

  virtual void AssembleMatrix(Matrix& A, Vector& u, NLInfo& nlinfo);
  virtual void AssembleMatrix(Matrix& A, Vector& u);
  virtual void ComputeIlu(Matrix& A, Vector& u);

  virtual void BoundaryInit(Vector& u) const;

  virtual void SolutionTransfer(Vector& u) const;
  virtual void Transfer(Vector& u) const;

  virtual void vmulteq(const Matrix& A, Vector& y, const Vector& x) const;

  virtual void LinearMg(IndexType minlevel,
                        IndexType maxlevel,
                        const Matrix& A,
                        Vector& u,
                        const Vector& f,
                        CGInfo&);

  virtual double ComputeFunctional(Vector& f,
                                   const Vector& u,
                                   const std::string& label);

  virtual void AssembleDualMatrix(Matrix& A, Vector& u);

  // fuer gmres

  virtual void precondition(const Matrix& A, Vector& x, Vector& y);
  virtual void DeleteVector(Vector& p);
  virtual void Equ(Vector& dst, double s, const Vector& src) const;
  virtual void Zero(Vector& dst) const;

  virtual void AddNodeVector(const std::string& name, Vector& q);
  virtual void DeleteNodeVector(const std::string& q);

  virtual void newton(Matrix& A,
                      Vector& u,
                      const Vector& f,
                      Vector& r,
                      Vector& w,
                      NLInfo& info);
};
} // namespace Gascoigne

#endif
