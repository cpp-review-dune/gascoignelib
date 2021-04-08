/**
 *
 * Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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

#ifndef __StdMultiLevelSolverData_h
#define __StdMultiLevelSolverData_h

#include "cginfo.h"
#include "paramfile.h"
#include <iostream>
#include <map>

/**********************************************************/

namespace Gascoigne {
class StdMultiLevelSolverData {
protected:
  std::string _solver, _mgtype, _linearsolve, _nonlinearsolve;
  mutable int _countresidual;
  int _coarselevel;
  int _i_show_nonlinear_comp_residuals, _i_show_linear_comp_residuals,
      _i_show_comp_residual_names;
  int _i_save_nonlinear_comp_residuals, _i_save_linear_comp_residuals;
  double _mgomega;

  int _gmresmemsize;
  CGInfo precinfo;

public:
  StdMultiLevelSolverData() {}
  StdMultiLevelSolverData(const ParamFile &param);
  virtual ~StdMultiLevelSolverData();

  virtual void BasicInit(const ParamFile &param);

  std::string &Solver() { return _solver; }
  int &CountResidual() const { return _countresidual; }

  int CoarseLevel() const { return _coarselevel; }

  int &CoarseLevel() { return _coarselevel; }

  std::string &MgType() { return _mgtype; }
  double &MgOmega() { return _mgomega; }
  const std::string &LinearSolve() const { return _linearsolve; }
  std::string &NonLinearSolve() { return _nonlinearsolve; }

  int GmresMemSize() const { return _gmresmemsize; }
  CGInfo &GetPrecInfo() { return precinfo; }
  const CGInfo &GetPrecInfo() const { return precinfo; }

  int SaveNonLinearCompResiduals() const {
    return _i_save_nonlinear_comp_residuals;
  }
  int SaveLinearCompResiduals() const { return _i_save_linear_comp_residuals; }
  int ShowNonLinearCompResiduals() const {
    return _i_show_nonlinear_comp_residuals;
  }
  int ShowLinearCompResiduals() const { return _i_show_linear_comp_residuals; }
  int ShowCompResidualNames() const {
    return (_i_show_comp_residual_names && (_i_show_nonlinear_comp_residuals ||
                                            _i_show_linear_comp_residuals));
  }
};

/**********************************************************/

} // namespace Gascoigne

/**********************************************************/

#endif
