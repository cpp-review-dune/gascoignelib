/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#ifndef  __NavierStokesGls3d_h
#define  __NavierStokesGls3d_h

#include  "navierstokes3d.h"
#include  "glsequation.h"
#include  "glsstabilization.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class NavierStokesGls3d : public NavierStokes3d, public virtual GlsEquation
{
  protected:

  mutable GlsStabilization ST;
  
  public:

  ~NavierStokesGls3d();
  NavierStokesGls3d();
  NavierStokesGls3d(const ParamFile* pf);

  std::string GetName() const;

  void SetTime(double time, double dt) const {Application::SetTime(time,dt); ST.DeltaT() = dt;}
  //
  /// Computation of gls stabilization parameters
  //
  void glspoint(double h, const FemFunction& U, const Vertex3d& v) const;

  //
  /// for Galerkin-Least-Squares
  //
  void L(DoubleVector& dst, const FemFunction& U) const;
  void S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;

  void LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;
};
}

#endif
