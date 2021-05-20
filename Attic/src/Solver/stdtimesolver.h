/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008 by the Gascoigne 3D authors
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

#ifndef __StdTimeSolver_h
#define __StdTimeSolver_h

#include "stdsolver.h"

/*-------------------------------------------------------------*/

namespace Gascoigne {

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear Solver for time-dependent Equations

///
///
//////////////////////////////////////////////

class StdTimeSolver : public virtual StdSolver
{
private:
  TimePattern _TP;
  MatrixInterface* _MMP;

protected:
  double _dt, _theta, _time;
  std::array<double, 2> _rhs;

  MatrixInterface*& GetMassMatrixPointer() { return _MMP; }

  const TimePattern& GetTimePattern() const { return _TP; }
  TimePattern& GetTimePattern() { return _TP; }

  virtual MatrixInterface* NewMassMatrix(int ncomp,
                                         const std::string& matrixtype);
  virtual std::string PrecondCGMass(GlobalVector& u,
                                    GlobalVector& f,
                                    const TimePattern& TP,
                                    double s = 1.);

public:
  StdTimeSolver();
  ~StdTimeSolver();

  void RegisterMatrix();
  void ReInitMatrix();

  void SetTimeData(double dt,
                   double theta,
                   double time,
                   double oldrhs = -1.,
                   double newrhs = 1.);
  double GetTime() const { return _time; }
  virtual GascoigneVisualization* NewGascoigneVisualization() const;

  void SetProblem(const ProblemDescriptorInterface& PDX);

  void InitialCondition(Vector& f, double d = 1.) const;
  void TimeRhsOperator(Vector& f, const Vector& u) const;
  void TimeRhs(int k, Vector& f) const;
  void Form(Vector& y, const Vector& x, double d) const;
  void AssembleMatrix(Matrix& A, const Vector& u, double d);
  std::string GetName() const;
  void L2Projection(Vector& u, Vector& f);

  void SetMassMatrix(MatrixInterface& MM, bool init = false);
  const MatrixInterface* GetMassMatrix() const { return _MMP; }
  MatrixInterface* GetMassMatrix() { return _MMP; }
  void MassMatrixVector(Vector& gf, const Vector& gu, double d) const;
};
} // namespace Gascoigne

#endif
