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

#ifndef __StdLoop_h
#define __StdLoop_h

#include "adaptordata.h"
#include "basicloop.h"
#include "extrapolator.h"

/*-----------------------------------------*/

namespace Gascoigne {
//////////////////////////////////////////////
//
///@brief
/// Implementation of diverse outer loops: mesh adaptation, time stepping...

///
///
//////////////////////////////////////////////

class StdLoop : public virtual BasicLoop
{
protected:
  // for computational time measurement
  mutable StopWatch _clock_newmesh, _clock_solve, _clock_write;

  int _nmin, _nmax, _coarse;
  double _p;
  int _random_coarsening;

  /// if yes: print statistic on the runtime after each iteration
  bool _runtime_statistics;

  std::string _estimator, _extrapolate, _refiner;
  DoubleVector _JErr;
  Extrapolator Extra;

  // new vectors

  DoubleVector ComputeFunctionals(Vector& f, Vector& u);

  const DoubleVector GetExactValues() const;
  const std::vector<std::string> GetFunctionalNames() const;

  virtual void EtaVisu(std::string name, int i, const DoubleVector& eta) const;
  virtual void EtaCellVisu(std::string name,
                           int i,
                           const GlobalVector& eta) const;
  virtual void AdaptMesh(const DoubleVector& eta);
  virtual void AdaptMesh(const DoubleVector& eta,
                         std::string refine_or_coarsen_step);
  virtual DoubleVector Functionals(Vector& u, Vector& f, bool output = true);
  virtual double Estimator(DoubleVector& eta, Vector& u, Vector& f);

public:
  StdLoop();
  StdLoop(const ParamFile& paramfile,
          const ProblemContainer* PC,
          const FunctionalContainer* FC);
  ~StdLoop();

  void BasicInit(const ParamFile& paramfile,
                 const ProblemContainer* PC,
                 const FunctionalContainer* FC);

  void run(const std::string& problemlabel);
  void ClockOutput() const;
};
} // namespace Gascoigne
/*-----------------------------------------*/

#endif
