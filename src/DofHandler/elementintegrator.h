/*----------------------------   elementintegrator.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef __elementintegrator_H
#define __elementintegrator_H
/*----------------------------   elementintegrator.h
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

#include "basicintegrator.h"
#include "integrationformula.h"
#include "integrationformulasummed.h"
#include "triaintegrationformula.h"
#include "patchintegrationformula.h"

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  ... comments ElementIntegrator
////
////
/////////////////////////////////////////////

template<int DIM, class IFF, class IFE, class IFB, class IFM>
class ElementIntegrator : public BasicIntegrator
{
private:
protected:
  double Volume2MeshSize(double vol) const
  {
    return pow(vol, 1. / static_cast<double>(DIM));
  }

public:
  //
  ////  Con(De)structor
  //
  ElementIntegrator(){};
  ~ElementIntegrator(){};

  // Integration-Formulas
  virtual IntegrationFormulaInterface* NewIFF() const { return new IFF; }
  virtual IntegrationFormulaInterface* NewIFE() const { return new IFE; }
  virtual IntegrationFormulaInterface* NewIFB() const { return new IFB; }
  virtual IntegrationFormulaInterface* NewIFM() const { return new IFM; }

  std::string GetName() const { return "ElementIntegrator"; }

  void Rhs(const DomainRightHandSide& RHS,
           LocalVector& F,
           const FemInterface& FEM,
           const LocalData& Q,
           const LocalData& QC) const;
  void Form(const Equation& EQ,
            LocalVector& F,
            const FemInterface& FEM,
            const LocalVector& U,
            const LocalData& Q,
            const LocalData& QC) const;
  void AdjointForm(const Equation& EQ,
                   LocalVector& F,
                   const FemInterface& FEM,
                   const LocalVector& U,
                   const LocalData& Q,
                   const LocalData& QC) const;
  void BoundaryForm(const BoundaryEquation& BE,
                    LocalVector& F,
                    const FemInterface& FEM,
                    const LocalVector& U,
                    int ile,
                    int col,
                    const LocalData& Q,
                    const LocalData& QC) const;
  void Matrix(const Equation& EQ,
              EntryMatrix& E,
              const FemInterface& FEM,
              const LocalVector& U,
              const LocalData& Q,
              const LocalData& QC) const;
  void BoundaryMatrix(const BoundaryEquation& BE,
                      EntryMatrix& E,
                      const FemInterface& FEM,
                      const LocalVector& U,
                      int ile,
                      int col,
                      const LocalData& Q,
                      const LocalData& QC) const;
  double MassMatrix(EntryMatrix& E, const FemInterface& FEM) const;
  void BoundaryMassMatrix(EntryMatrix& E,
                          const FemInterface& FEM,
                          int ile) const;
  void MassForm(const TimePattern& TP,
                LocalVector& F,
                const FemInterface& FEM,
                const LocalVector& U) const;

  void RhsPoint(LocalVector& b,
                const FemInterface& E,
                const Vertex<DIM>& p,
                int comp) const;
  void DiracRhsPoint(LocalVector& b,
                     const FemInterface& E,
                     const Vertex<DIM>& p,
                     const DiracRightHandSide& DRHS,
                     int j,
                     const LocalData& Q,
                     const LocalData& QC) const;

  ////////////////////////////////////////////////// Functionals
  double LocalDiv(const FemInterface& F, const LocalVector& U) const;

  double ComputePointValue(const FemInterface& E,
                           const Vertex<DIM>& p,
                           const LocalVector& U,
                           int comp) const;
  double ComputeDomainFunctional(const DomainFunctional& F,
                                 const FemInterface& FEM,
                                 const LocalVector& U,
                                 const LocalData& Q,
                                 const LocalData& QC) const;
  double ComputeErrorDomainFunctional(const DomainFunctional& F,
                                      const FemInterface& FEM,
                                      const LocalVector& U,
                                      const LocalData& Q,
                                      const LocalData& QC) const;
  double ComputeBoundaryFunctional(const BoundaryFunctional& F,
                                   const FemInterface& FEM,
                                   int ile,
                                   int col,
                                   const LocalVector& U,
                                   const LocalData& Q,
                                   const LocalData& QC) const;
  void EvaluateCellRightHandSide(LocalVector& F,
                                 const DomainRightHandSide& CF,
                                 const FemInterface& FEM,
                                 const LocalData& Q,
                                 const LocalData& QC) const;
  void EvaluateBoundaryCellRightHandSide(LocalVector& F,
                                         const BoundaryRightHandSide& CF,
                                         const FemInterface& FEM,
                                         int ile,
                                         int col,
                                         const LocalData& Q,
                                         const LocalData& QC) const;

  void ErrorsByExactSolution(LocalVector& dst,
                             const FemInterface& FE,
                             const ExactSolution& ES,
                             const LocalVector& U,
                             const LocalData& Q,
                             const LocalData& QC) const;

  void BoundaryRhs(const BoundaryRightHandSide& RHS,
                   LocalVector& F,
                   const FemInterface& FEM,
                   int ile,
                   int col,
                   const LocalData& Q,
                   const LocalData& QC) const;

  void IntegrateMassDiag(DoubleVector& F, const FemInterface& FEM) const;

  void IntegrateBoundaryMassDiag(DoubleVector& F,
                                 const FemInterface& FEM,
                                 int ile,
                                 int col) const;

  // no warning for overloaded virtual function (2d/3d)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
  void RhsCurve(LocalVector& F,
                const FemInterface& FEM,
                Vertex<DIM>& xr0,
                Vertex<DIM>& xr1,
                double H,
                double ND0,
                double ND1,
                int ncomp,
                int comp) const;
#pragma GCC diagnostic pop
};

#define ElementIntegratorP12d                                                  \
  ElementIntegrator<2,                                                         \
                    TriaQuadFormula<TriaSixpointFormula>,                      \
                    TriaQuadFormula<TriaSixpointFormula>,                      \
                    LineGauss2,                                                \
                    TriaQuadFormula<TriaSixpointFormula>>

#define ElementIntegratorQ12d                                                  \
  ElementIntegrator<2, QuadGauss4, QuadGauss9, LineGauss2, QuadGauss4>
#define ElementIntegratorQ22d                                                  \
  ElementIntegrator<2, QuadGauss9, QuadGauss16, LineGauss3, QuadGauss9>
#define ElementIntegratorQ42d                                                  \
  ElementIntegrator<2, QuadGauss25, QuadGauss36, LineGauss5, QuadGauss25>

#define ElementIntegratorQ13d                                                  \
  ElementIntegrator<3, HexGauss8, HexGauss27, QuadGauss4, HexGauss8>
#define ElementIntegratorQ23d                                                  \
  ElementIntegrator<3, HexGauss27, HexGauss64, QuadGauss9, HexGauss27>
#define ElementIntegratorQ43d                                                  \
  ElementIntegrator<3, HexGauss125, HexGauss216, QuadGauss16, HexGauss125>

typedef ElementIntegrator<2,
			  PatchFormula2d<4, QuadGauss4>,
			  PatchFormula2d<9, QuadGauss9>,
			  PatchFormula1d<2, LineGauss2>,
			  PatchFormula2d<4, QuadGauss4>>
  ElementIntegratorQ12dPatch;


  
} // namespace Gascoigne

/*----------------------------   elementintegrator.h
 * ---------------------------*/
/* end of #ifndef __elementintegrator_H */
#endif
/*----------------------------   elementintegrator.h
 * ---------------------------*/
