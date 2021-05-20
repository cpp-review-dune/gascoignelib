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

#ifndef __MultiLevelAlgorithm_h
#define __MultiLevelAlgorithm_h

#include "algorithm.h"
#include "multilevelsolver.h"

/*-----------------------------------------*/

namespace Gascoigne {

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class MultiLevelAlgorithm : public Algorithm
{
private:
  MultiLevelSolver* _S;

  int _coarselevel;

  std::string _mgtype, _linearsolve;
  double _mgomega;

protected:
  virtual const SolverInterface* GetSolver() const { return _S->GetSolver(); }
  virtual SolverInterface* GetSolver() { return _S->GetSolver(); }

  virtual const SolverInterface* GetSolver(int i) const
  {
    return _S->GetSolver(i);
  }
  virtual SolverInterface* GetSolver(int i) { return _S->GetSolver(i); }

  virtual const MultiLevelSolver* GetMultiLevelSolver() const { return _S; }
  virtual MultiLevelSolver* GetMultiLevelSolver() { return _S; }

  void DeleteVector(Vector& u) const { _S->DeleteVector(u); }
  void AssembleMatrixAndIlu(Vector& u);
  void LinearSolve(Vector& du, const Vector& y, CGInfo& cginfo);
  void NonLinear(Vector& u,
                 Vector& f,
                 const std::string& problemlabel,
                 int iter);
  void VWCycle(std::vector<double>& res,
               std::vector<double>& rw,
               int l,
               int finelevel,
               int coarselevel,
               const std::string& p,
               Vector& u,
               Vector& b,
               Vector& v);
  void LinearMGSolve(Vector& du, const Vector& y, CGInfo& cginfo);
  virtual void Precondition(Vector& x, Vector& y);

public:
  MultiLevelAlgorithm()
    : _S(NULL)
  {}
  virtual ~MultiLevelAlgorithm() {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
  void BasicInit(const ParamFile* paramfile,
                 MultiLevelSolver* MLS,
                 const NumericInterface* NI,
                 const ProblemContainer* PC);
#pragma GCC diagnostic pop

  void ReInitVector(Vector& u) const { _S->ReInitVector(u); }
  void RunLinear(const std::string& problemlabel);
  void RunNonLinear(const std::string& problemlabel, int iter = 0);
  void GlobalRefineLoop(const std::string& problemlabel);
  void LocalRefineLoop(const std::string& problemlabel,
                       FunctionalContainer* FC = NULL);
};
} // namespace Gascoigne

/*-----------------------------------------*/

#endif
