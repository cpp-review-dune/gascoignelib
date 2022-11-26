/*----------------------------   elementintegrator.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef __mixedelementintegrator_H
#define __mixedelementintegrator_H
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

/**
 * The mixed element integrator uses a different finite element
 * for the test- and the trial-functions.
 * Both hoewever do rely on the same dof-handling.
 * The typical option:
 *   BaseQ1Patch for the trial
 *   BaseQ2      for the test-space
 **/

#include "../Discretization/Q1/basicintegrator.h"
#include "../Discretization/Q1/integrationformula.h"
#include "../Discretization/Q1/integrationformulasummed.h"
#include "../Discretization/Q1/triaintegrationformula.h"
#include "../Discretization/Q2/patchintegrationformula.h"

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  ... comments ElementIntegrator
////
////
/////////////////////////////////////////////

template<int DIM, class IFF, class IFE, class IFB, class IFM>
class MixedElementIntegrator : public BasicIntegrator
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
  MixedElementIntegrator(){};
  ~MixedElementIntegrator(){};

  // Integration-Formulas
  virtual IntegrationFormulaInterface* NewIFF() const { return new IFF; }
  virtual IntegrationFormulaInterface* NewIFE() const { return new IFE; }
  virtual IntegrationFormulaInterface* NewIFB() const { return new IFB; }
  virtual IntegrationFormulaInterface* NewIFM() const { return new IFM; }

  std::string GetName() const { return "MixedElementIntegrator"; }

  void Rhs(const DomainRightHandSide& RHS,
           LocalVector& F,
           const FemInterface& FEMTRIAL,
           const FemInterface& FEMTEST,
           const LocalData& Q,
           const LocalData& QC) const;
  void Form(const Equation& EQ,
            LocalVector& F,
            const FemInterface& FEMTRIAL,
            const FemInterface& FEMTEST,
            const LocalVector& U,
            const LocalData& Q,
            const LocalData& QC) const;
  void AdjointForm(const Equation& EQ,
                   LocalVector& F,
                   const FemInterface& FEM,
                   const LocalVector& U,
                   const LocalData& Q,
                   const LocalData& QC) const
  {
    abort();
  }
  void BoundaryForm(const BoundaryEquation& BE,
                    LocalVector& F,
                    const FemInterface& FEM,
                    const LocalVector& U,
                    int ile,
                    int col,
                    const LocalData& Q,
                    const LocalData& QC) const
  {
    abort();
  }
  void Matrix(const Equation& EQ,
              EntryMatrix& E,
              const FemInterface& FEMTRIAL,
              const FemInterface& FEMTEST,
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
                      const LocalData& QC) const
  {
    abort();
  }
  double MassMatrix(EntryMatrix& E, const FemInterface& FEM) const { abort(); }
  void BoundaryMassMatrix(EntryMatrix& E,
                          const FemInterface& FEM,
                          int ile) const
  {
    abort();
  }
  void MassForm(const TimePattern& TP,
                LocalVector& F,
                const FemInterface& FEM,
                const LocalVector& U) const
  {
    abort();
  }

  void RhsPoint(LocalVector& b,
                const FemInterface& E,
                const Vertex<DIM>& p,
                int comp) const
  {
    abort();
  }
  void DiracRhsPoint(LocalVector& b,
                     const FemInterface& E,
                     const Vertex<DIM>& p,
                     const DiracRightHandSide& DRHS,
                     int j,
                     const LocalData& Q,
                     const LocalData& QC) const
  {
    abort();
  }

  ////////////////////////////////////////////////// Functionals
  double LocalDiv(const FemInterface& F, const LocalVector& U) const
  {
    abort();
  }

  double ComputePointValue(const FemInterface& E,
                           const Vertex<DIM>& p,
                           const LocalVector& U,
                           int comp) const
  {
    abort();
  }
  double ComputeDomainFunctional(const DomainFunctional& F,
                                 const FemInterface& FEM,
                                 const LocalVector& U,
                                 const LocalData& Q,
                                 const LocalData& QC) const
  {
    abort();
  }
  double ComputeErrorDomainFunctional(const DomainFunctional& F,
                                      const FemInterface& FEM,
                                      const LocalVector& U,
                                      const LocalData& Q,
                                      const LocalData& QC) const
  {
    abort();
  }
  double ComputeBoundaryFunctional(const BoundaryFunctional& F,
                                   const FemInterface& FEM,
                                   int ile,
                                   int col,
                                   const LocalVector& U,
                                   const LocalData& Q,
                                   const LocalData& QC) const
  {
    abort();
  }
  void EvaluateCellRightHandSide(LocalVector& F,
                                 const DomainRightHandSide& CF,
                                 const FemInterface& FEM,
                                 const LocalData& Q,
                                 const LocalData& QC) const
  {
    abort();
  }
  void EvaluateBoundaryCellRightHandSide(LocalVector& F,
                                         const BoundaryRightHandSide& CF,
                                         const FemInterface& FEM,
                                         int ile,
                                         int col,
                                         const LocalData& Q,
                                         const LocalData& QC) const
  {
    abort();
  }

  void ErrorsByExactSolution(LocalVector& dst,
                             const FemInterface& FE,
                             const ExactSolution& ES,
                             const LocalVector& U,
                             const LocalData& Q,
                             const LocalData& QC) const
  {
    abort();
  }

  void BoundaryRhs(const BoundaryRightHandSide& RHS,
                   LocalVector& F,
                   const FemInterface& FEM,
                   int ile,
                   int col,
                   const LocalData& Q,
                   const LocalData& QC) const
  {
    abort();
  }

  void IntegrateMassDiag(DoubleVector& F, const FemInterface& FEM) const
  {
    abort();
  }

  void IntegrateBoundaryMassDiag(DoubleVector& F,
                                 const FemInterface& FEM,
                                 int ile,
                                 int col) const
  {
    abort();
  }

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
                int comp) const
  {
    abort();
  }
#pragma GCC diagnostic pop
};

#define MixedElementIntegratorQ12dPatch                                        \
  MixedElementIntegrator<2,                                                    \
                         PatchFormula2d<9, QuadGauss9>,                        \
                         PatchFormula2d<16, QuadGauss16>,                      \
                         PatchFormula1d<3, LineGauss3>,                        \
                         PatchFormula2d<9, QuadGauss9>>

} // namespace Gascoigne

/*----------------------------   elementintegrator.h
 * ---------------------------*/
/* end of #ifndef __elementintegrator_H */
#endif
/*----------------------------   elementintegrator.h
 * ---------------------------*/
