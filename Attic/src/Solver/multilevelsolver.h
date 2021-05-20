/**
 *
 * Copyright (C) 2008, 2009 by the Gascoigne 3D authors
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

#ifndef __MultiLevelSolver_h
#define __MultiLevelSolver_h

#include "functionalcontainer.h"
#include "meshagentinterface.h"
#include "mginterpolatorinterface.h"
#include "numericinterface.h"
#include "paramfile.h"
#include "problemcontainer.h"
#include "problemdescriptorinterface.h"
#include "solverinterface.h"

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

class MultiLevelSolver
{
protected:
  std::vector<SolverInterface*> _SP;
  std::vector<MgInterpolatorInterface*> _Interpolator;

  const MeshAgentInterface* _MAP;
  const ParamFile* _paramfile;
  const ProblemContainer* _PC;
  const ProblemDescriptorInterface* _PD;
  const NumericInterface* _NI;

  int ComputeLevel;

  const MeshAgentInterface* GetMeshAgent() const { return _MAP; }

  void NewMgInterpolator();
  void SolverNewMesh();
  void ReInitMatrix();
  void RegisterMatrix();
  virtual void NewSolvers();
  // virtual SolverInterface* NewSolver(int solverlevel);
  void Transfer(int high, int low, Vector& u) const;
  void SolutionTransfer(Vector& u) const;

public:
  // Constructor

  MultiLevelSolver();
  virtual ~MultiLevelSolver();

  std::string GetName() const { return "MultiLevelSolver"; }
  int nlevels() const
  {
    assert(GetMeshAgent());
    return GetMeshAgent()->nlevels();
  }

  const ProblemContainer* GetProblemContainer() const { return _PC; };

  int FinestLevel() const { return nlevels() - 1; }

  SolverInterface* GetSolver(int l)
  {
    assert(l < _SP.size());
    return _SP[l];
  }
  const SolverInterface* GetSolver(int l) const
  {
    assert(l < _SP.size());
    return _SP[l];
  }
  SolverInterface* GetSolver()
  {
    assert(_SP.size() == nlevels());
    return _SP[FinestLevel()];
  }
  const SolverInterface* GetSolver() const
  {
    assert(_SP.size() == nlevels());
    return _SP[FinestLevel()];
  }

  std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers()
  {
    return _Interpolator;
  }
  const std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers() const
  {
    return _Interpolator;
  }

  void BasicInit(const NumericInterface* NI,
                 const MeshAgentInterface*,
                 const ParamFile* paramfile,
                 const ProblemContainer* PC);

  void ReInitVector(Vector& v);
  void DeleteVector(Vector& g);
  void SetProblem(const std::string& label);
  void ReInit(const std::string& problemlabel);
  void AssembleMatrix(Vector& u);
  void AssembleDualMatrix(Vector& u);
  void ComputeIlu();
  void ComputeIlu(Vector& u);
  const DoubleVector ComputeFunctionals(Vector& f,
                                        const Vector& u,
                                        FunctionalContainer* FC) const;
};
} // namespace Gascoigne

#endif
