/*----------------------------   navierstokes.h ---------------------------*/
/*      $Id:$                 */
#ifndef __navierstokes_H
#define __navierstokes_H
/*----------------------------   navierstokes.h ---------------------------*/

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
 * \brief Data class for the Navier-Stokes equations, stationary and
 * non-stationary
 *
 * Provides the required persistent data to fully specify the Navier-Stokes
 * equations
 *
 * visc:            viscosity (always mandatory)
 * symmetrictensor: 1 (nabla v + nabla v^T), 0 (nabla v)
 *
 * alpha: Parameter for the LPS pressure stabilization
 * delta: Parameter for the LPS convection stabilization
 *
 * theta:  Parameter for the theta time-stepping scheme
 * dt,time: time step and current time
 *
 * outflowcolor: color of the outflow do-nothing boundary
 */
class NavierStokesData {
public:
  double visc, alpha, delta;

  bool symmetrictensor;
  int outflowcolor;
  double dt, theta;

  // reads parameters from the input file
  void BasicInit(const ParamFile &pf) {
    DataFormatHandler DFH;
    DFH.insert("visc", &visc, 1.);   // viscosity, must be positive
    DFH.insert("alpha", &alpha, 0.); // lps pressure stabilization
    DFH.insert("delta", &delta, 0.); // lps convection stabilization
    DFH.insert("theta", &theta, 0.); // theta - time stepping
    DFH.insert("dt", &dt, 0.);       // time step
    DFH.insert("symmetrictensor", &symmetrictensor, 0); // use full tensor?
    DFH.insert("outflowcolor", &outflowcolor, 0); // do-nothing boundary color

    FileScanner FS(DFH, pf, "Equation");
    assert(visc > 0);
  }
};

////////////////////////////////////////////////// Stationary

/**
 *
 * \brief Navier-Stokes equations, 2d and 3d, LPS stabilization
 *
 */
template <int DIM>
class NavierStokesLps : public virtual LpsEquation, public BoundaryEquation {

protected:
  /// Data object
  NavierStokesData data;

  // current normal vector at the boundary
  mutable Vertex<DIM> normal;
  // local lps parameters
  mutable double alpha, delta;

public:
  NavierStokesLps<DIM>(const NavierStokesData &PD) : data(PD) {}

  NavierStokesLps<DIM> *createNew() const {
    return new NavierStokesLps<DIM>(data);
  }

  std::string GetName() const { return "NavierStokes - LPS - "; }

  int GetNcomp() const { return DIM + 1; }

  /// Variational formulation
  void Form(VectorIterator b, const FemFunction &U,
            const TestFunction &N) const {
    for (int i = 0; i < DIM; ++i) {
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
        b[i + 1] += U[j + 1].m() * U[i + 1][j + 1] * N.m();
    }
  }

  /// Jacobian
  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N) const {
    for (int i = 0; i < DIM; ++i) {
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

      // convection
      for (int j = 0; j < DIM; ++j) {
        A(i + 1, i + 1) += U[j + 1].m() * M[j + 1] * N.m();
        A(i + 1, j + 1) += M.m() * U[i + 1][j + 1] * N.m();
      }
    }
  }

  ////////////////////////////////////////////////// Correction of the the
  /// do-nothing boundary

  /// store the normal in each quadrature point on the boundary
  void pointboundary(double h, const FemFunction &U, const Vertex<DIM> &v,
                     const Vertex<DIM> &n) const {
    normal = n;
  }

  /// variational form on the boundary
  void Form(VectorIterator b, const FemFunction &U, const TestFunction &N,
            int col) const {
    if (data.symmetrictensor ==
        false) // only correct outflow boundary if symmetric tensor is used
      return;

    if (col == data.outflowcolor) {
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
          b[i + 1] -= data.visc * normal[j] * U[j + 1][i + 1] * N.m();
    }
  }

  /// Boundary-Jacobian
  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N, int col) const {
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
  void lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const {
    double vel = 0;

    alpha = data.alpha / (data.visc / h / h + vel / h);
    delta = data.delta / (data.visc / h / h + vel / h);
  }

  /// Variational form of LPS
  void StabForm(VectorIterator b, const FemFunction &U, const FemFunction &UP,
                const TestFunction &NP) const {
    // pressure stabilization
    for (int i = 0; i < DIM; ++i)
      b[0] += alpha * UP[0][i + 1] * NP[i + 1];

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
  void StabMatrix(EntryMatrix &A, const FemFunction &U, const TestFunction &Np,
                  const TestFunction &Mp) const {
    // pressure stabilization
    for (int i = 0; i < DIM; ++i)
      A(0, 0) += alpha * Mp[i + 1] * Np[i + 1];

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

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

/**
 *
 * \brief Navier-Stokes equations, 2d and 3d, LPS stabilization, Non-STationary
 *
 */
template <int DIM>
class NavierStokesLpsTime : public virtual LpsEquation,
                            public BoundaryEquation {

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
  mutable FemFunction *OLD;

public:
  NavierStokesLpsTime<DIM>(const NavierStokesData &PD) : data(PD) {
    assert(PD.theta > 0 && PD.theta <= 1.0); // theta parameter must be in (0,1]
  }

  NavierStokesLpsTime<DIM> *createNew() const {
    return new NavierStokesLpsTime<DIM>(data);
  }

  std::string GetName() const { return "NavierStokes - LPS - Time"; }

  int GetNcomp() const { return DIM + 1; }

  /// Retrieves old solution OLD
  void SetFemData(FemData &q) const {
    assert(q.find("OLD") != q.end());
    OLD = &q["OLD"];
  }

  /**
   * \brief Variational formulation
   *
   * We implement the theta time stepping scheme, taking the pressure
   * and the divergence equation fully implicit
   */
  void Form(VectorIterator b, const FemFunction &U,
            const TestFunction &N) const {
    for (int i = 0; i < DIM; ++i) {
      // Time derivative
      b[i + 1] += (U[i + 1].m() - (*OLD)[i + 1].m()) / data.dt * N.m();

      // divergence equation
      b[0] += U[i + 1][i + 1] * N.m();

      // viscosity
      for (int j = 0; j < DIM; ++j)
        b[i + 1] += data.theta * data.visc * U[i + 1][j + 1] * N[j + 1];

      if (data.symmetrictensor)
        for (int j = 0; j < DIM; ++j)
          b[i + 1] += data.theta * data.visc * U[j + 1][i + 1] * N[j + 1];

      // pressure
      b[i + 1] -= U[0].m() * N[i + 1];

      // convection
      for (int j = 0; j < DIM; ++j)
        b[i + 1] += data.theta * U[j + 1].m() * U[i + 1][j + 1] * N.m();

      // add explicit part
      if (data.theta < 1) {
        // viscosity
        for (int j = 0; j < DIM; ++j)
          b[i + 1] +=
              (1.0 - data.theta) * data.visc * (*OLD)[i + 1][j + 1] * N[j + 1];

        if (data.symmetrictensor)
          for (int j = 0; j < DIM; ++j)
            b[i + 1] += (1.0 - data.theta) * data.visc * (*OLD)[j + 1][i + 1] *
                        N[j + 1];

        // convection
        for (int j = 0; j < DIM; ++j)
          b[i + 1] += (1.0 - data.theta) * (*OLD)[j + 1].m() *
                      (*OLD)[i + 1][j + 1] * N.m();
      }
    }
  }

  /// Jacobian
  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N) const {
    for (int i = 0; i < DIM; ++i) {
      // Time derivative
      A(i + 1, i + 1) += M.m() / data.dt * N.m();

      // divergence equation
      A(0, i + 1) += M[i + 1] * N.m();

      // viscosity
      for (int j = 0; j < DIM; ++j)
        A(i + 1, i + 1) += data.theta * data.visc * M[j + 1] * N[j + 1];
      if (data.symmetrictensor)
        for (int j = 0; j < DIM; ++j)
          A(i + 1, j + 1) += data.theta * data.visc * M[i + 1] * N[j + 1];

      // pressure
      A(i + 1, 0) -= M.m() * N[i + 1];

      // convection
      for (int j = 0; j < DIM; ++j) {
        A(i + 1, i + 1) += data.theta * U[j + 1].m() * M[j + 1] * N.m();
        A(i + 1, j + 1) += data.theta * M.m() * U[i + 1][j + 1] * N.m();
      }
    }
  }

  ////////////////////////////////////////////////// Correction of the the
  /// do-nothing boundary

  /// store the normal in each quadrature point on the boundary
  void pointboundary(double h, const FemFunction &U, const Vertex<DIM> &v,
                     const Vertex<DIM> &n) const {
    normal = n;
  }

  /// variational form on the boundary
  void Form(VectorIterator b, const FemFunction &U, const TestFunction &N,
            int col) const {
    if (data.symmetrictensor ==
        false) // only correct outflow boundary if symmetric tensor is used
      return;

    if (col == data.outflowcolor) {
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j) {
          b[i + 1] -=
              data.theta * data.visc * normal[j] * U[j + 1][i + 1] * N.m();
          b[i + 1] -= (1.0 - data.theta) * data.visc * normal[j] *
                      (*OLD)[j + 1][i + 1] * N.m();
        }
    }
  }

  /// Boundary-Jacobian
  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N, int col) const {
    if (data.symmetrictensor ==
        false) // only correct outflow boundary if symmetric tensor is used
      return;

    if (col == data.outflowcolor) {
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
          A(i + 1, j + 1) -=
              data.theta * data.visc * normal[j] * M[i + 1] * N.m();
    }
  }

  ////////////////////////////////////////////////// LPS Stabilization

  /// define local stabilization parameters
  void lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const {
    double vel = 0;

    alpha = data.alpha / (data.visc / h / h + vel / h + 1.0 / data.dt);
    delta = data.delta / (data.visc / h / h + vel / h + 1.0 / data.dt);
  }

  /// Variational form of LPS
  void StabForm(VectorIterator b, const FemFunction &U, const FemFunction &UP,
                const TestFunction &NP) const {
    // pressure stabilization
    for (int i = 0; i < DIM; ++i)
      b[0] += alpha * UP[0][i + 1] * NP[i + 1];

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
  void StabMatrix(EntryMatrix &A, const FemFunction &U, const TestFunction &Np,
                  const TestFunction &Mp) const {
    // pressure stabilization
    for (int i = 0; i < DIM; ++i)
      A(0, 0) += alpha * Mp[i + 1] * Np[i + 1];

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

#define NavierStokesLps2d NavierStokesLps<2>
#define NavierStokesLps3d NavierStokesLps<3>

#define NavierStokesLpsTime2d NavierStokesLpsTime<2>
#define NavierStokesLpsTime3d NavierStokesLpsTime<3>

} // namespace Gascoigne

/*----------------------------   navierstokes.h ---------------------------*/
/* end of #ifndef __navierstokes_H */
#endif
/*----------------------------   navierstokes.h ---------------------------*/
