/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 by the Gascoigne 3D authors
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

#ifndef __DiscretizationInterface_h
#define __DiscretizationInterface_h

#include "boundaryequation.h"
#include "boundaryfunctional.h"
#include "boundaryrighthandside.h"
#include "compvector.h"
#include "curve.h"
#include "datacontainer.h"
#include "diracrighthandside.h"
#include "dirichletdata.h"
#include "domainfunction.h"
#include "domainfunctional.h"
#include "domainrighthandside.h"
#include "equation.h"
#include "exactsolution.h"
#include "gascoigne.h"
#include "gascoignemesh.h"
#include "matrixinterface.h"
#include "meshtransferinterface.h"
#include "mginterpolatorinterface.h"
#include "paramfile.h"
#include "pointfunctional.h"
#include "problemdescriptorbase.h"
#include <string>

namespace Gascoigne {

/////////////////////////////////////////////
///
///@brief
///  ... comments DiscretizationInterface

///
///
/////////////////////////////////////////////

class DiscretizationInterface
{
private:
protected:
public:
  DiscretizationInterface() {}
  virtual ~DiscretizationInterface() {}

  virtual const DataContainer& GetDataContainer() const = 0;
  virtual void SetDataContainer(const DataContainer& q) = 0;

  //
  //// Functions called from the Solver
  //
  virtual std::string GetName() const = 0;

  virtual void AddNodeVector(const std::string& name,
                             const GlobalVector* q) const = 0;
  virtual void DeleteNodeVector(const std::string& name) const = 0;

  virtual void AddCellVector(const std::string& name,
                             const GlobalVector* q) const = 0;
  virtual void DeleteCellVector(const std::string& name) const = 0;

  virtual void AddParameterVector(const std::string& name,
                                  const GlobalParameterVector* q) const = 0;
  virtual void DeleteParameterVector(const std::string& name) const = 0;

  virtual void BasicInit(const ParamFile& pf) = 0;
  virtual void ReInit(const GascoigneMesh* M) = 0;

  virtual int ndofs() const = 0; // returns the number of degrees of freedom
  virtual int nelements()
    const = 0; // returns the number of elements in the discretization
  virtual int ndegree() const
  {
    std::cerr << "\"DiscretizationInterface::ndegree\" not written!\n"
                 "Continuing anyways."
              << std::endl;
    return 1;
  }
  virtual int ndofs_withouthanging() const { return ndofs(); }
  virtual Vertex2d vertex2d(int i) const { abort(); }
  virtual Vertex3d vertex3d(int i) const { abort(); }

  virtual void Structure(SparseStructureInterface* S) const = 0;
  virtual void Form(GlobalVector& f,
                    const GlobalVector& u,
                    const Equation& EQ,
                    double d) const
  {
    assert(0);
  }
  virtual void Rhs(GlobalVector& f,
                   const DomainRightHandSide& RHS,
                   double s) const
  {
    assert(0);
  };
  virtual void Matrix(MatrixInterface& A,
                      const GlobalVector& u,
                      const Equation& EQ,
                      double) const
  {
    assert(0);
  }

  virtual void AdjointForm(GlobalVector& f,
                           const GlobalVector& u,
                           const Equation& EQ,
                           double d) const
  {
    std::cerr << "\"DiscretizationInterface::AdjointForm\" not written!"
              << std::endl;
    abort();
  }
  virtual void BoundaryMatrix(MatrixInterface& A,
                              const GlobalVector& u,
                              const IntSet& Colors,
                              const BoundaryEquation& BE,
                              double d) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryMatrix\" not written!"
              << std::endl;
    abort();
  }

  // Visualization
  virtual void VisuVtk(const ComponentInformation* CI,
                       const ParamFile& pf,
                       const std::string& name,
                       const GlobalVector& u,
                       int i) const
  {
    std::cerr << "\"DiscretizationInterface::VisuVtk not written!" << std::endl;
    abort();
  }

  // New Inteface.
  virtual void BoundaryForm(GlobalVector& f,
                            const GlobalVector& u,
                            const ProblemDescriptorInterface& PD,
                            double d) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryForm-NEW\" not written!"
              << std::endl;
    abort();
  }
  virtual void BoundaryMatrix(MatrixInterface& A,
                              const GlobalVector& u,
                              const ProblemDescriptorInterface& PD,
                              double d) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryMatrix-NEW\" not written!"
              << std::endl;
    abort();
  }

  virtual void MassMatrix(MatrixInterface& M) const
  {
    std::cerr << "\"DiscretizationInterface::MassMatrix\" not written!"
              << std::endl;
    abort();
  }

  virtual void BoundaryMassMatrix(MatrixInterface& A,
                                  const IntSet& Colors) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryMassMatrix\" not written!"
              << std::endl;
    abort();
  }
  virtual void MassForm(GlobalVector& f,
                        const GlobalVector& u,
                        const TimePattern& TP,
                        double s) const
  {
    std::cerr << "\"DiscretizationInterface::MassForm\" not written!"
              << std::endl;
    abort();
  }
  virtual void DiracRhs(GlobalVector& f,
                        const DiracRightHandSide& DRHS,
                        double s) const
  {
    std::cerr << "\"DiscretizationInterface::DiracRhs\" not written!"
              << std::endl;
    abort();
  }
  virtual void BoundaryRhs(GlobalVector& f,
                           const IntSet& Colors,
                           const BoundaryRightHandSide& BRHS,
                           double s) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryRhs\" not written!"
              << std::endl;
    abort();
  }
  virtual void HNAverage(GlobalVector& x) const {}
  virtual void HNDistribute(GlobalVector& x) const {}
  virtual void HNZero(GlobalVector& x) const {}
  virtual bool HNZeroCheck(const GlobalVector& x) const { return false; }
  virtual void HNAverageData() const {}
  virtual void HNZeroData() const {}
  virtual void Interpolate(GlobalVector& u,
                           const DomainInitialCondition& U) const
  {
    std::cerr << "\"DiscretizationInterface::Interpolate\" not written!"
              << std::endl;
    abort();
  }
  virtual void InterpolateSolution(GlobalVector& u,
                                   const GlobalVector& uold) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateSolution\" not written!"
              << std::endl;
    abort();
  }
  virtual void InterpolateDirac(GlobalVector& u, const GlobalVector& uold) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateDirac\" not written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletMatrix(MatrixInterface& A,
                                     int col,
                                     const std::vector<int>& comp) const
  {
    std::cerr
      << "\"DiscretizationInterface::StrongDirichletmatrix\" not written!"
      << std::endl;
    abort();
  }
  virtual void StrongDirichletMatrixOnlyRow(MatrixInterface& A,
                                            int col,
                                            const std::vector<int>& comp) const
  {
    std::cerr << "\"DiscretizationInterface::StrongDirichletMatrixOnlyRow\" "
                 "not written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletVector(GlobalVector& u,
                                     const DirichletData& BF,
                                     int col,
                                     const std::vector<int>& comp,
                                     double d = 1.) const
  {
    std::cerr
      << "\"DiscretizationInterface::StronDirichletVector\" not written!"
      << std::endl;
    abort();
  }
  virtual void StrongDirichletVectorZero(GlobalVector& u,
                                         int col,
                                         const std::vector<int>& comp) const
  {
    std::cerr
      << "\"DiscretizationInterface::StrongDirichletVectorZero\" not written!"
      << std::endl;
    abort();
  }
  virtual void StrongPeriodicVector(GlobalVector& u,
                                    const PeriodicData& BF,
                                    int col,
                                    const std::vector<int>& comp,
                                    double d = 1.) const
  {
    std::cerr
      << "\"DiscretizationInterface::StrongPeriodicVector\" not written!"
      << std::endl;
    abort();
  }

  ////////////////////////////////////////////////// New Interface Colors and
  /// Components in ProblemDescriptor
  virtual void StrongDirichletVector(GlobalVector& u,
                                     const DirichletData* DD,
                                     double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::StronDirichletVector - NEW\" not "
                 "written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletVectorZero(
    GlobalVector& u,
    const ProblemDescriptorInterface& PD) const
  {
    std::cerr << "\"DiscretizationInterface::StrongDirichletVectorZero - NEW\" "
                 "not written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletMatrix(MatrixInterface& A,
                                     const ProblemDescriptorInterface& PD) const
  {
    std::cerr << "\"DiscretizationInterface::StrongDirichletmatrix - NEW\" not "
                 "written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletMatrixOnlyRow(
    MatrixInterface& A,
    const ProblemDescriptorInterface& PD) const
  {
    std::cerr << "\"DiscretizationInterface::StrongDirichletMatrixOnlyRow - "
                 "NEW\" not written!"
              << std::endl;
    abort();
  }

  virtual void InitFilter(DoubleVector&) const
  {
    std::cerr << "\"DiscretizationInterface::InitFilter\" not written!"
              << std::endl;
    abort();
  }

  virtual void StabForm(GlobalVector& f,
                        const GlobalVector& u,
                        const ProblemDescriptorInterface& PD,
                        double d) const
  {
    std::cerr << "\"DiscretizationInterface::StabForm\" not written!"
              << std::endl;
    abort();
  }

  // Functionals
  virtual void ComputeError(const GlobalVector& u,
                            LocalVector& err,
                            const ExactSolution* ES) const
  {
    std::cerr << "\"DiscretizationInterface::ComputeError\" not written!"
              << std::endl;
    abort();
  }
  virtual void AssembleError(GlobalVector& eta,
                             const GlobalVector& u,
                             LocalVector& err,
                             const ExactSolution* ES) const
  {
    std::cerr << "\"DiscretizationInterface::AssembleError\" not written!"
              << std::endl;
    abort();
  }

  ////////////////////////////////////////////////// Functionals
  virtual double LocalDiv(const LocalVector& U, const LocalVector& M) const
  {
    return 0.0;
  }

  virtual double ComputeBoundaryFunctional(const GlobalVector& u,
                                           const IntSet& Colors,
                                           const BoundaryFunctional& BF) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputeBoundaryFunctional\" not written!"
      << std::endl;
    abort();
  }
  virtual double ComputeDomainFunctional(const GlobalVector& u,
                                         const DomainFunctional& F) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputeDomainFunctional\" not written!"
      << std::endl;
    abort();
  }
  virtual double ComputeErrorDomainFunctional(const GlobalVector& u,
                                              const DomainFunctional& F) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputeDomainFunctional\" not written!"
      << std::endl;
    abort();
  }

  virtual double ComputePointFunctional(const GlobalVector& u,
                                        const PointFunctional& FP) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputePointFunctional\" not written!"
      << std::endl;
    abort();
  }

  virtual void EvaluateCellRightHandSide(GlobalVector& f,
                                         const DomainRightHandSide& CF,
                                         double d = 1.) const
  {
    std::cerr
      << "\"DiscretizationInterface::EvaluateCellRighthandside\" not written!"
      << std::endl;
    abort();
  }

  virtual void EvaluateBoundaryCellRightHandSide(
    GlobalVector& f,
    const IntSet& Colors,
    const BoundaryRightHandSide& CF,
    double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::"
                 "EvaluateBoundaryCellRighthandside\" not written!"
              << std::endl;
    abort();
  }

  virtual void EvaluateParameterRightHandSide(GlobalVector& f,
                                              const DomainRightHandSide& CF,
                                              double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::EvaluateParameterRighthandside\" "
                 "not written!"
              << std::endl;
    abort();
  }

  virtual void EvaluateBoundaryParameterRightHandSide(
    GlobalVector& f,
    const IntSet& Colors,
    const BoundaryRightHandSide& CF,
    double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::"
                 "EvaluateBoundaryParameterRighthandside\" not written!"
              << std::endl;
    abort();
  }

  virtual void InterpolateDomainFunction(GlobalVector& f,
                                         const DomainFunction& DF) const
  {
    std::cerr
      << "\"DiscretizationInterface::InterpolateDomainFunction\" not written!"
      << std::endl;
    abort();
  }

  virtual void InterpolateCellDomainFunction(GlobalVector& f,
                                             const DomainFunction& DF) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateCellDomainFunction\" "
                 "not written!"
              << std::endl;
    abort();
  }

  virtual void ConstructInterpolator(MgInterpolatorInterface* I,
                                     const MeshTransferInterface* MT)
  {
    std::cerr
      << "\"DiscretizationInterface::ConstructInterpolator\" not written!"
      << std::endl;
    abort();
  }

  virtual void GetVolumes(DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetVolumes\" not written!"
              << std::endl;
    abort();
  }

  virtual void GetAreas(DoubleVector& a, const IntSet& Colors) const
  {
    std::cerr << "\"DiscretizationInterface::GetAreas\" not written!"
              << std::endl;
    abort();
  }

  virtual void GetMassDiag(DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetMassDiag\" not written!"
              << std::endl;
    abort();
  }

  virtual void GetBoundaryMassDiag(DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetBoundaryMassDiag\" not written!"
              << std::endl;
    abort();
  }

  virtual void RhsCurve(GlobalVector& F, const Curve& C, int comp, int N) const
  {
    std::cerr << "\"DiscretizationInterface::RhsCurve\" not written!"
              << std::endl;
    abort();
  }
};
} // namespace Gascoigne

#endif
