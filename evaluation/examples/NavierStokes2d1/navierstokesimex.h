/*----------------------------   navierstokesimex.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __navierstokesimex_H
#define __navierstokesimex_H
/*----------------------------   navierstokesimex.h     ---------------------------*/


/**
 *
 * Copyright (C) 2004,2018,2020 by the Gascoigne 3D authors
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

#include "boundaryequation.h"
#include "equation.h"
#include "filescanner.h"
#include "lpsequation.h"
#include "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne {

/**
 *
 * \brief Navier-Stokes equations, 2d and 3d, LPS stabilization, Non-STationary
 *
 */
template<int DIM>
class NavierStokesLpsTimeIMEX
  : public virtual LpsEquation
  , public BoundaryEquation
{

protected:
  /// Data object
  NavierStokesData data;

  // current normal vector at the boundary
  mutable Vertex<DIM> normal;
  // local lps parameters
  mutable double alpha, delta;
  /**
   * Solution at last time step. This solution is submitted in the Loop using
   * GetMultiLevelSolver()->AddNodeVector("OLD",old)
   * Here, the keyword OLD is found in the function SetFemData
   */
  mutable FemFunction *OLD, *OLDOLD;

public:
  NavierStokesLpsTimeIMEX<DIM>(const NavierStokesData& PD)
    : data(PD)
  {
  }

  NavierStokesLpsTimeIMEX<DIM>* createNew() const
  {
    return new NavierStokesLpsTimeIMEX<DIM>(data);
  }

  std::string GetName() const { return "NavierStokes - LPS - Time IMEX"; }

  int GetNcomp() const { return DIM + 1; }

  /// Retrieves old solution OLD
  void SetFemData(FemData& q) const
  {
    assert(q.find("OLD") != q.end());
    OLD = &q["OLD"];
    assert(q.find("OLDOLD") != q.end());
    OLDOLD = &q["OLDOLD"];
  }

  /**
   * \brief Variational formulation
   *
   * We implement the BDF2 time stepping scheme, taking the pressure
   * and the divergence equation fully implicit.
   * The convection term is extrapolated via (2 OLD - OLDOLD)
   */
  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    for (int i = 0; i < DIM; ++i) {
      // Time derivative
      b[i + 1] += 0.5/data.dt * (3.0 * U[i + 1].m() - 4.0 * (*OLD)[i + 1].m() + (*OLDOLD)[i+1].m()) * N.m();

      // divergence equation
      b[0] += U[i + 1][i + 1] * N.m();

      // viscosity
      for (int j = 0; j < DIM; ++j)
        b[i + 1] += data.visc * U[i + 1][j + 1] * N[j + 1];

      if (data.symmetrictensor)
        for (int j = 0; j < DIM; ++j)
          b[i + 1] += data.visc * U[j + 1][i + 1] * N[j + 1];

      // pressure
      b[i + 1] -= U[0].m() * N[i + 1];

      // convection
      for (int j = 0; j < DIM; ++j)
        b[i + 1] += (2.0 * (*OLD)[j + 1].m() - (*OLDOLD)[j+1].m()) *
	  (2.0 * (*OLD)[i + 1][j + 1] - (*OLDOLD)[i+1][j+1]) * N.m();
    }
  }

  /// Jacobian
  void Matrix(EntryMatrix& A,
              const FemFunction& U,
              const TestFunction& M,
              const TestFunction& N) const
  {
    for (int i = 0; i < DIM; ++i) {
      // Time derivative
      A(i + 1, i + 1) += 1.5 * M.m() / data.dt * N.m();

      // divergence equation
      A(0, i + 1) += M[i + 1] * N.m();

      // viscosity
      for (int j = 0; j < DIM; ++j)
        A(i + 1, i + 1) += data.visc * M[j + 1] * N[j + 1];
      if (data.symmetrictensor)
        for (int j = 0; j < DIM; ++j)
          A(i + 1, j + 1) += data.visc * M[i + 1] * N[j + 1];

      // pressure
      A(i + 1, 0) -= M.m() * N[i + 1];
    }
  }

  ////////////////////////////////////////////////// Correction of the the
  /// do-nothing boundary

  /// store the normal in each quadrature point on the boundary
  void pointboundary(double h,
                     const FemFunction& U,
                     const Vertex<DIM>& v,
                     const Vertex<DIM>& n) const
  {
    normal = n;
  }

  /// variational form on the boundary
  void Form(VectorIterator b,
            const FemFunction& U,
            const TestFunction& N,
            int col) const
  {
    if (data.symmetrictensor ==
        false) // only correct outflow boundary if symmetric tensor is used
      return;
    
    if (col == data.outflowcolor) {
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j) {
          b[i + 1] -=data.visc * normal[j] * U[j + 1][i + 1] * N.m();
        }
    }
  }
  
  /// Boundary-Jacobian
  void Matrix(EntryMatrix& A,
              const FemFunction& U,
              const TestFunction& M,
              const TestFunction& N,
              int col) const
  {
    if (data.symmetrictensor ==
        false) // only correct outflow boundary if symmetric tensor is used
      return;

    if (col == data.outflowcolor) {
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
          A(i + 1, j + 1) -= data.visc * normal[j] * M[i + 1] * N.m();
    }
  }

  ////////////////////////////////////////////////// LPS Stabilization

  /// define local stabilization parameters
  void lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    double vel = 0;

    alpha = data.alpha / (data.visc / h / h + vel / h + 1.0 / data.dt);
    delta = data.delta / (data.visc / h / h + vel / h + 1.0 / data.dt);
  }

  /// Variational form of LPS
  void StabForm(VectorIterator b,
                const FemFunction& U,
                const FemFunction& UP,
                const TestFunction& NP) const
  {
    // pressure stabilization
    for (int i = 0; i < DIM; ++i)
      b[0] += alpha * UP[0][i + 1] * NP[i + 1];

    return;
    
    
    // convective stabilization
    if (data.delta > 0) {
      double convNP = 0;
      for (int j = 0; j < DIM; ++j)
        convNP += U[j + 1].m() * NP[j + 1];
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
          b[i + 1] += delta * U[j + 1].m() * UP[i + 1][j + 1] * convNP;
    }
  }

  /// Jacobian of LPS
  void StabMatrix(EntryMatrix& A,
                  const FemFunction& U,
                  const TestFunction& Np,
                  const TestFunction& Mp) const
  {
    // pressure stabilization
    for (int i = 0; i < DIM; ++i)
      A(0, 0) += alpha * Mp[i + 1] * Np[i + 1];

    return;
    
    // convective stabilization. Only derivatives w.r.t. UP are computed
    if (data.delta > 0) {
      double convNP = 0;
      for (int j = 0; j < DIM; ++j)
        convNP += U[j + 1].m() * Np[j + 1];
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
          A(i + 1, i + 1) += delta * U[j + 1].m() * Mp[j + 1] * convNP;
    }
  }
};

} // namespace Gascoigne


/*----------------------------   navierstokesimex.h     ---------------------------*/
/* end of #ifndef __navierstokesimex_H */
#endif
/*----------------------------   navierstokesimex.h     ---------------------------*/
