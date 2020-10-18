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
#include "mginterpolatorinterface.h"
#include "monitor.h"
#include "problemcontainer.h"
#include "problemdescriptorinterface.h"
#include "stdmultilevelsolverdata.h"
#include "stopwatch.h"
#include "meshagentinterface.h"
#include "stdsolver.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
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
protected:  ///!!!!
  std::vector<StdSolver*> _SP;
  const MeshAgentInterface* _MAP;
  std::vector<MgInterpolatorInterface*> _Interpolator;

  mutable StopWatch _clock_residual, _clock_solve;
  mutable int ComputeLevel;
  mutable int oldnlevels;

  ParamFile _paramfile;

  Monitor MON;
  StdMultiLevelSolverData DataP;
  const ProblemDescriptorInterface* _PD;
  const ProblemContainer* _PC;
  const FunctionalContainer* _FC;

protected:
public:
  //////////////////////////////////////////////////
  // Constructor, Init

  StdMultiLevelSolver();
  explicit StdMultiLevelSolver(const MeshAgentInterface* GMGM, const ParamFile& paramfile,
                      const ProblemContainer* PC, const FunctionalContainer* FC = NULL);
  virtual ~StdMultiLevelSolver();

  virtual std::string GetName() const
  {
    return "StdMultiLevelSolver";
  }

  virtual void BasicInit(const MeshAgentInterface* GMGM, const ParamFile& paramfile,
                         const ProblemContainer* PC,
                         const FunctionalContainer* FC = NULL);

  virtual void ReInitMatrix(const Matrix& A);
  virtual void ReInitVector(VectorInterface& v);
  virtual void ReInitVector(VectorInterface& v, int comp);

  //////////////////////////////////////////////////
  // Access

  virtual const MeshAgentInterface* GetMeshAgent() const
  {
    return _MAP;
  }
  virtual std::vector<StdSolver*>& GetSolverPointers()
  {
    return _SP;
  }
  virtual const std::vector<StdSolver*>& GetSolverPointers() const
  {
    return _SP;
  }
  virtual const ProblemDescriptorInterface* GetProblemDescriptor() const
  {
    return _PD;
  }
  virtual StdSolver*& GetSolverPointer(int l)
  {
    assert(l < _SP.size());
    return _SP[l];
  }

  //////////////////////////////////////////////////
  // Solver
  virtual void NewSolvers();

  virtual StdSolver* NewSolver(int solverlevel);
  virtual void NewMgInterpolator();
  virtual void SolverNewMesh();
  virtual void SetComputeLevel(int level)
  {
    ComputeLevel = level;
  }
  virtual std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers()
  {
    return _Interpolator;
  }
  virtual const std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers() const
  {
    return _Interpolator;
  }
  virtual const ProblemContainer* GetProblemContainer() const
  {
    assert(_PC);
    return _PC;
  }
  virtual void SetProblemContainer(const ProblemContainer* PC)
  {
    _PC = PC;
  }
  virtual const FunctionalContainer* GetFunctionalContainer() const
  {
    return _FC;
  }
  virtual void SetFunctionalContainer(const FunctionalContainer* FC)
  {
    _FC = FC;
  }

  virtual StdSolver* GetSolver(int l)
  {
    assert(l < _SP.size());
    return _SP[l];
  }
  virtual const StdSolver* GetSolver(int l) const
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
  virtual const DoubleVector ComputeFunctionals(VectorInterface& f,
                                                const VectorInterface& u);
  virtual const DoubleVector ComputeFunctionals(VectorInterface& f,
                                                const VectorInterface& u,
                                                FunctionalContainer* FC);

  //////////////////////////////////////////////////
  // Newton and Multigrid
  virtual double NewtonNorm(const VectorInterface& u) const
  {
    return GetSolver(ComputeLevel)->NewtonNorm(u);
  }
  virtual void mgstep(std::vector<double>& res, std::vector<double>& rw, int l, int maxl,
                      int minl, std::string& p0, std::string p,
		      const Matrix& A, VectorInterface& u,
                      VectorInterface& b, VectorInterface& v);

  virtual void Cg(VectorInterface& x, const VectorInterface& f, CGInfo& info);
  virtual void Gmres(const Matrix& A, VectorInterface& x, const VectorInterface& f, CGInfo& info);

  virtual void SolutionTransfer(int high, int low, VectorInterface& u) const;
  virtual void Transfer(int high, int low, VectorInterface& u) const;

  virtual void ViewProtocoll() const;

  
  virtual int nlevels() const
  {
    assert(GetMeshAgent());
    return GetMeshAgent()->nlevels();
  }
  virtual int FinestLevel() const
  {
    return nlevels() - 1;
  }
  virtual int CoarsestLevel() const
  {
    return 0;
  }

  Monitor& GetMonitor()
  {
    return MON;
  }

  const Monitor& GetMonitor() const
  {
    return MON;
  }

  // ReInit is to be called, whenever the mesh (hierarchy) has changed
  virtual void ReInit();
  // SetProblem announces the current ProblemDescriptor to all solvers in the hierarchy
  virtual void SetProblem(const std::string& label);

  // neue vektoren

  virtual std::string LinearSolve(int level,
				  const Matrix& A, VectorInterface& u, const VectorInterface& b,
                                  CGInfo& info);
  virtual std::string Solve(int level,
			    Matrix& A, VectorInterface& x, const VectorInterface& b,
                            NLInfo& nlinfo);
  virtual std::string Solve(Matrix& A, VectorInterface& x, const VectorInterface& b, NLInfo& nlinfo)
  {
    return Solve(nlevels() - 1, A, x, b, nlinfo);
  }

  virtual void InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const;
  virtual void InterpolateCellSolution(VectorInterface& u,
                                       const GlobalVector& uold) const;

  virtual void NewtonVectorZero(VectorInterface& w) const;
  virtual double NewtonResidual(VectorInterface& y, const VectorInterface& x,
                                const VectorInterface& b) const;
  virtual double NewtonUpdate(double& rr, VectorInterface& x, VectorInterface& dx,
                              VectorInterface& r, const VectorInterface& f,
                              NLInfo& nlinfo);
  virtual void NewtonLinearSolve(const Matrix& A, VectorInterface& x, const VectorInterface& b,
                                 CGInfo& info);
  virtual void NewtonMatrixControl(Matrix& A, VectorInterface& u, NLInfo& nlinfo);
  virtual void NewtonOutput(NLInfo& nlinfo) const;
  virtual void NewtonPreProcess(VectorInterface& u, const VectorInterface& f,
                                NLInfo& info) const;
  virtual void NewtonPostProcess(VectorInterface& u, const VectorInterface& f,
                                 NLInfo& info) const;

  virtual void AssembleMatrix(Matrix& A, VectorInterface& u, NLInfo& nlinfo);
  virtual void AssembleMatrix(Matrix& A, VectorInterface& u);
  virtual void ComputeIlu(Matrix& A, VectorInterface& u);

  virtual void BoundaryInit(VectorInterface& u) const;

  virtual void SolutionTransfer(VectorInterface& u) const;
  virtual void Transfer(VectorInterface& u) const;

  virtual void vmulteq(const Matrix& A,
		       VectorInterface& y, const VectorInterface& x) const;

  virtual void LinearMg(int minlevel, int maxlevel,
			const Matrix& A, VectorInterface& u,
                        const VectorInterface& f, CGInfo&);

  virtual double ComputeFunctional(VectorInterface& f, const VectorInterface& u,
                                   const std::string& label);

  virtual void AssembleDualMatrix(Matrix& A, VectorInterface& u);

  // fuer gmres

  virtual void precondition(const Matrix& A, VectorInterface& x, VectorInterface& y);
  virtual void DeleteVector(VectorInterface& p);
  virtual void Equ(VectorInterface& dst, double s, const VectorInterface& src) const;
  virtual void Zero(VectorInterface& dst) const;

  virtual void AddNodeVector(const std::string& name, VectorInterface& q);
  virtual void DeleteNodeVector(const std::string& q);

  virtual void newton(Matrix& A,
		      VectorInterface& u, const VectorInterface& f, VectorInterface& r,
                      VectorInterface& w, NLInfo& info);
};
}  // namespace Gascoigne

#endif
