/*----------------------------   elementlpsintegrator.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef __elementlpsintegrator_H
#define __elementlpsintegrator_H
/*----------------------------   elementlpsintegrator.h
 * ---------------------------*/

/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2018 by the Gascoigne 3D authors
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

#include "elementintegrator.h"
#include "lpsequation.h"

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  ... comments ElementIntegrator
////
////
/////////////////////////////////////////////

template<int DIM, class IFF, class IFE, class IFB, class IFM>
class ElementLpsIntegrator : public ElementIntegrator<DIM, IFF, IFE, IFB, IFM>
{
private:
protected:
public:
  //
  ////  Con(De)structor
  //

  ElementLpsIntegrator(){};
  ~ElementLpsIntegrator(){};

  std::string GetName() const { return "ElementLpsIntegrator"; }

  void Projection(FemFunction& NLPS, const FemInterface& FEM) const
  {
    for (int ii = 0; ii < FEM.n(); ii++) {
      FEM.init_test_functions(NLPS[ii], 1., ii);
    }
    if (DIM == 2) {
      NLPS[0].equ(-0.25, NLPS[4]);
      NLPS[2].equ(-0.25, NLPS[4]);
      NLPS[6].equ(-0.25, NLPS[4]);
      NLPS[8].equ(-0.25, NLPS[4]);

      NLPS[0].add(-0.5, NLPS[1], -0.5, NLPS[3]);
      NLPS[2].add(-0.5, NLPS[1], -0.5, NLPS[5]);
      NLPS[6].add(-0.5, NLPS[3], -0.5, NLPS[7]);
      NLPS[8].add(-0.5, NLPS[5], -0.5, NLPS[7]);
    } else if (DIM == 3) {
      NLPS[0].equ(-0.125, NLPS[13]);
      NLPS[2].equ(-0.125, NLPS[13]);
      NLPS[6].equ(-0.125, NLPS[13]);
      NLPS[8].equ(-0.125, NLPS[13]);
      NLPS[18].equ(-0.125, NLPS[13]);
      NLPS[20].equ(-0.125, NLPS[13]);
      NLPS[24].equ(-0.125, NLPS[13]);
      NLPS[26].equ(-0.125, NLPS[13]);

      NLPS[0].add(-0.5, NLPS[1], -0.5, NLPS[3], -0.5, NLPS[9]);
      NLPS[2].add(-0.5, NLPS[1], -0.5, NLPS[5], -0.5, NLPS[11]);
      NLPS[6].add(-0.5, NLPS[3], -0.5, NLPS[7], -0.5, NLPS[15]);
      NLPS[8].add(-0.5, NLPS[5], -0.5, NLPS[7], -0.5, NLPS[17]);
      NLPS[18].add(-0.5, NLPS[19], -0.5, NLPS[21], -0.5, NLPS[9]);
      NLPS[20].add(-0.5, NLPS[19], -0.5, NLPS[23], -0.5, NLPS[11]);
      NLPS[24].add(-0.5, NLPS[21], -0.5, NLPS[25], -0.5, NLPS[15]);
      NLPS[26].add(-0.5, NLPS[23], -0.5, NLPS[25], -0.5, NLPS[17]);

      NLPS[0].add(-0.25, NLPS[4], -0.25, NLPS[10], -0.25, NLPS[12]);
      NLPS[2].add(-0.25, NLPS[4], -0.25, NLPS[10], -0.25, NLPS[14]);
      NLPS[6].add(-0.25, NLPS[4], -0.25, NLPS[16], -0.25, NLPS[12]);
      NLPS[8].add(-0.25, NLPS[4], -0.25, NLPS[16], -0.25, NLPS[14]);
      NLPS[18].add(-0.25, NLPS[22], -0.25, NLPS[10], -0.25, NLPS[12]);
      NLPS[20].add(-0.25, NLPS[22], -0.25, NLPS[10], -0.25, NLPS[14]);
      NLPS[24].add(-0.25, NLPS[22], -0.25, NLPS[16], -0.25, NLPS[12]);
      NLPS[26].add(-0.25, NLPS[22], -0.25, NLPS[16], -0.25, NLPS[14]);
    }
  }

  void Form(const Equation& EQ,
            LocalVector& F,
            const FemInterface& FEM,
            const LocalVector& U,
            const LocalData& Q,
            const LocalData& QC) const
  {
    ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::Form(EQ, F, FEM, U, Q, QC);

    FemFunction NLPS(FEM.n());
    FemFunction MLPS(FEM.n());
    FemFunction UHP;

    const LpsEquation& LEQ = dynamic_cast<const LpsEquation&>(EQ);

    IFF IF;
    Vertex<DIM> x, xi;

    for (int k = 0; k < IF.n(); k++) {
      IF.xi(xi, k);
      FEM.point(xi);

      double vol = FEM.J();
      double weight = IF.w(k) * vol;
      double h =
        ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::Volume2MeshSize(vol);

      FEM.x(x);
      BasicIntegrator::universal_point(
        FEM, ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_UH, U);
      BasicIntegrator::universal_point(
        FEM, ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_QH, Q);

      LEQ.SetFemData(ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_QH);
      LEQ.lpspoint(h, ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_UH, x);

      Projection(NLPS, FEM); // fuellt NLPS
      BasicIntegrator::universal_point(UHP, U, NLPS);
      for (int i = 0; i < FEM.n(); i++)
        MLPS[i].equ(weight, NLPS[i]);

      for (int i = 0; i < FEM.n(); i++) {
        LEQ.StabForm(F.start(i),
                     ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_UH,
                     UHP,
                     MLPS[i]);
      }
    }
  }

  void Matrix(const Equation& EQ,
              EntryMatrix& E,
              const FemInterface& FEM,
              const LocalVector& U,
              const LocalData& Q,
              const LocalData& QC) const
  {
    ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::Matrix(EQ, E, FEM, U, Q, QC);

    assert(E.Ndof() == FEM.n());
    assert(E.Mdof() == FEM.n());
    assert(E.Ncomp() == U.ncomp());

    FemFunction NLPS(FEM.n());
    FemFunction MLPS(FEM.n());

    const LpsEquation& LEQ = dynamic_cast<const LpsEquation&>(EQ);

    IFF IF;

    Vertex<DIM> x, xi;
    for (int k = 0; k < IF.n(); k++) {
      IF.xi(xi, k);
      FEM.point(xi);

      double vol = FEM.J();
      double weight = IF.w(k) * vol;
      FEM.x(x);

      BasicIntegrator::universal_point(
        FEM, ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_UH, U);
      BasicIntegrator::universal_point(
        FEM, ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_QH, Q);
      double h =
        ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::Volume2MeshSize(vol);
      LEQ.SetFemData(ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_QH);
      LEQ.lpspointmatrix(h, ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_UH, x);

      Projection(NLPS, FEM);

      for (int i = 0; i < FEM.n(); i++)
        MLPS[i].equ(weight, NLPS[i]);

      for (int j = 0; j < FEM.n(); j++) {
        for (int i = 0; i < FEM.n(); i++) {
          E.SetDofIndex(i, j);
          LEQ.StabMatrix(E,
                         ElementIntegrator<DIM, IFF, IFE, IFB, IFM>::_UH,
                         NLPS[i],
                         MLPS[j]);
        }
      }
    }
  }
};
} // namespace Gascoigne

#define ElementLpsIntegratorQ12d                                               \
  ElementLpsIntegrator<2,                                                      \
                       PatchFormula2d<4, QuadGauss4>,                          \
                       PatchFormula2d<9, QuadGauss9>,                          \
                       PatchFormula1d<2, LineGauss2>,                          \
                       PatchFormula2d<4, QuadGauss4>>
#define ElementLpsIntegratorQ22d                                               \
  ElementLpsIntegrator<2, QuadGauss9, QuadGauss16, LineGauss3, QuadGauss9>
#define ElementLpsIntegratorQ13d                                               \
  ElementLpsIntegrator<3,                                                      \
                       PatchFormula3d<8, HexGauss8>,                           \
                       PatchFormula3d<27, HexGauss27>,                         \
                       PatchFormula2d<4, QuadGauss4>,                          \
                       PatchFormula3d<8, HexGauss8>>
#define ElementLpsIntegratorQ23d                                               \
  ElementLpsIntegrator<3, HexGauss27, HexGauss64, QuadGauss9, HexGauss27>

/*----------------------------   elementlpsintegrator.h
 * ---------------------------*/
/* end of #ifndef __elementlpsintegrator_H */
#endif
/*----------------------------   elementlpsintegrator.h
 * ---------------------------*/
