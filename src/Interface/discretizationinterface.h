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
  virtual Vertex2d vertex2d([[maybe_unused]] int i) const { abort(); }
  virtual Vertex3d vertex3d([[maybe_unused]] int i) const { abort(); }

  virtual void Structure(SparseStructureInterface* S) const = 0;
  virtual void Form([[maybe_unused]] GlobalVector& f,
                    [[maybe_unused]] const GlobalVector& u,
                    [[maybe_unused]] const Equation& EQ,
                    [[maybe_unused]] double d) const
  {
    assert(0);
  }
  virtual void Rhs([[maybe_unused]] GlobalVector& f,
                   [[maybe_unused]] const DomainRightHandSide& RHS,
                   [[maybe_unused]] double s) const
  {
    assert(0);
  };
  virtual void Matrix([[maybe_unused]] MatrixInterface& A,
                      [[maybe_unused]] const GlobalVector& u,
                      [[maybe_unused]] const Equation& EQ,
                      [[maybe_unused]] double) const
  {
    assert(0);
  }

  virtual void AdjointForm([[maybe_unused]] GlobalVector& f,
                           [[maybe_unused]] const GlobalVector& u,
                           [[maybe_unused]] const Equation& EQ,
                           [[maybe_unused]] double d) const
  {
    std::cerr << "\"DiscretizationInterface::AdjointForm\" not written!"
              << std::endl;
    abort();
  }
  virtual void BoundaryMatrix([[maybe_unused]] MatrixInterface& A,
                              [[maybe_unused]] const GlobalVector& u,
                              [[maybe_unused]] const IntSet& Colors,
                              [[maybe_unused]] const BoundaryEquation& BE,
                              [[maybe_unused]] double d) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryMatrix\" not written!"
              << std::endl;
    abort();
  }

  // Visualization
  virtual void VisuVtk([[maybe_unused]] const ComponentInformation* CI,
                       [[maybe_unused]] const ParamFile& pf,
                       [[maybe_unused]] const std::string& name,
                       [[maybe_unused]] const GlobalVector& u,
                       [[maybe_unused]] int i) const
  {
    std::cerr << "\"DiscretizationInterface::VisuVtk not written!" << std::endl;
    abort();
  }

  // New Inteface.
  virtual void BoundaryForm(
    [[maybe_unused]] GlobalVector& f,
    [[maybe_unused]] const GlobalVector& u,
    [[maybe_unused]] const ProblemDescriptorInterface& PD,
    [[maybe_unused]] double d) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryForm-NEW\" not written!"
              << std::endl;
    abort();
  }
  virtual void BoundaryMatrix(
    [[maybe_unused]] MatrixInterface& A,
    [[maybe_unused]] const GlobalVector& u,
    [[maybe_unused]] const ProblemDescriptorInterface& PD,
    [[maybe_unused]] double d) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryMatrix-NEW\" not written!"
              << std::endl;
    abort();
  }

  virtual void MassMatrix([[maybe_unused]] MatrixInterface& M) const
  {
    std::cerr << "\"DiscretizationInterface::MassMatrix\" not written!"
              << std::endl;
    abort();
  }

  virtual void BoundaryMassMatrix([[maybe_unused]] MatrixInterface& A,
                                  [[maybe_unused]] const IntSet& Colors) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryMassMatrix\" not written!"
              << std::endl;
    abort();
  }
  virtual void MassForm([[maybe_unused]] GlobalVector& f,
                        [[maybe_unused]] const GlobalVector& u,
                        [[maybe_unused]] const TimePattern& TP,
                        [[maybe_unused]] double s) const
  {
    std::cerr << "\"DiscretizationInterface::MassForm\" not written!"
              << std::endl;
    abort();
  }
  virtual void DiracRhs([[maybe_unused]] GlobalVector& f,
                        [[maybe_unused]] const DiracRightHandSide& DRHS,
                        [[maybe_unused]] double s) const
  {
    std::cerr << "\"DiscretizationInterface::DiracRhs\" not written!"
              << std::endl;
    abort();
  }
  virtual void BoundaryRhs([[maybe_unused]] GlobalVector& f,
                           [[maybe_unused]] const IntSet& Colors,
                           [[maybe_unused]] const BoundaryRightHandSide& BRHS,
                           [[maybe_unused]] double s) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryRhs\" not written!"
              << std::endl;
    abort();
  }
  virtual void HNAverage([[maybe_unused]] GlobalVector& x) const {}
  virtual void HNDistribute([[maybe_unused]] GlobalVector& x) const {}
  virtual void HNZero([[maybe_unused]] GlobalVector& x) const {}
  virtual bool HNZeroCheck([[maybe_unused]] const GlobalVector& x) const
  {
    return false;
  }
  virtual void HNAverageData() const {}
  virtual void HNZeroData() const {}
  virtual void Interpolate(
    [[maybe_unused]] GlobalVector& u,
    [[maybe_unused]] const DomainInitialCondition& U) const
  {
    std::cerr << "\"DiscretizationInterface::Interpolate\" not written!"
              << std::endl;
    abort();
  }
  virtual void InterpolateSolution(
    [[maybe_unused]] GlobalVector& u,
    [[maybe_unused]] const GlobalVector& uold) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateSolution\" not written!"
              << std::endl;
    abort();
  }
  virtual void InterpolateDirac([[maybe_unused]] GlobalVector& u,
                                [[maybe_unused]] const GlobalVector& uold) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateDirac\" not written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletMatrix(
    [[maybe_unused]] MatrixInterface& A,
    [[maybe_unused]] int col,
    [[maybe_unused]] const std::vector<int>& comp) const
  {
    std::cerr
      << "\"DiscretizationInterface::StrongDirichletmatrix\" not written!"
      << std::endl;
    abort();
  }
  virtual void StrongDirichletMatrixOnlyRow(
    [[maybe_unused]] MatrixInterface& A,
    [[maybe_unused]] int col,
    [[maybe_unused]] const std::vector<int>& comp) const
  {
    std::cerr << "\"DiscretizationInterface::StrongDirichletMatrixOnlyRow\" "
                 "not written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletVector(
    [[maybe_unused]] GlobalVector& u,
    [[maybe_unused]] const DirichletData& BF,
    [[maybe_unused]] int col,
    [[maybe_unused]] const std::vector<int>& comp,
    [[maybe_unused]] double d = 1.) const
  {
    std::cerr
      << "\"DiscretizationInterface::StronDirichletVector\" not written!"
      << std::endl;
    abort();
  }
  virtual void StrongDirichletVectorZero(
    [[maybe_unused]] GlobalVector& u,
    [[maybe_unused]] int col,
    [[maybe_unused]] const std::vector<int>& comp) const
  {
    std::cerr
      << "\"DiscretizationInterface::StrongDirichletVectorZero\" not written!"
      << std::endl;
    abort();
  }
  virtual void StrongPeriodicVector(
    [[maybe_unused]] GlobalVector& u,
    [[maybe_unused]] const PeriodicData& BF,
    [[maybe_unused]] int col,
    [[maybe_unused]] const std::vector<int>& comp,
    [[maybe_unused]] double d = 1.) const
  {
    std::cerr
      << "\"DiscretizationInterface::StrongPeriodicVector\" not written!"
      << std::endl;
    abort();
  }

  ////////////////////////////////////////////////// New Interface Colors and
  /// Components in ProblemDescriptor
  virtual void StrongDirichletVector([[maybe_unused]] GlobalVector& u,
                                     [[maybe_unused]] const DirichletData* DD,
                                     [[maybe_unused]] double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::StronDirichletVector - NEW\" not "
                 "written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletVectorZero(
    [[maybe_unused]] GlobalVector& u,
    [[maybe_unused]] const ProblemDescriptorInterface& PD) const
  {
    std::cerr << "\"DiscretizationInterface::StrongDirichletVectorZero - NEW\" "
                 "not written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletMatrix(
    [[maybe_unused]] MatrixInterface& A,
    [[maybe_unused]] const ProblemDescriptorInterface& PD) const
  {
    std::cerr << "\"DiscretizationInterface::StrongDirichletmatrix - NEW\" not "
                 "written!"
              << std::endl;
    abort();
  }
  virtual void StrongDirichletMatrixOnlyRow(
    [[maybe_unused]] MatrixInterface& A,
    [[maybe_unused]] const ProblemDescriptorInterface& PD) const
  {
    std::cerr << "\"DiscretizationInterface::StrongDirichletMatrixOnlyRow - "
                 "NEW\" not written!"
              << std::endl;
    abort();
  }

  virtual void InitFilter([[maybe_unused]] DoubleVector&) const
  {
    std::cerr << "\"DiscretizationInterface::InitFilter\" not written!"
              << std::endl;
    abort();
  }

  virtual void StabForm([[maybe_unused]] GlobalVector& f,
                        [[maybe_unused]] const GlobalVector& u,
                        [[maybe_unused]] const ProblemDescriptorInterface& PD,
                        [[maybe_unused]] double d) const
  {
    std::cerr << "\"DiscretizationInterface::StabForm\" not written!"
              << std::endl;
    abort();
  }

  // Functionals
  virtual void ComputeError([[maybe_unused]] const GlobalVector& u,
                            [[maybe_unused]] LocalVector& err,
                            [[maybe_unused]] const ExactSolution* ES) const
  {
    std::cerr << "\"DiscretizationInterface::ComputeError\" not written!"
              << std::endl;
    abort();
  }
  virtual void AssembleError([[maybe_unused]] GlobalVector& eta,
                             [[maybe_unused]] const GlobalVector& u,
                             [[maybe_unused]] LocalVector& err,
                             [[maybe_unused]] const ExactSolution* ES) const
  {
    std::cerr << "\"DiscretizationInterface::AssembleError\" not written!"
              << std::endl;
    abort();
  }

  ////////////////////////////////////////////////// Functionals
  virtual double LocalDiv([[maybe_unused]] const LocalVector& U,
                          [[maybe_unused]] const LocalVector& M) const
  {
    return 0.0;
  }

  virtual double ComputeBoundaryFunctional(
    [[maybe_unused]] const GlobalVector& u,
    [[maybe_unused]] const IntSet& Colors,
    [[maybe_unused]] const BoundaryFunctional& BF) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputeBoundaryFunctional\" not written!"
      << std::endl;
    abort();
  }
  virtual double ComputeDomainFunctional(
    [[maybe_unused]] const GlobalVector& u,
    [[maybe_unused]] const DomainFunctional& F) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputeDomainFunctional\" not written!"
      << std::endl;
    abort();
  }
  virtual double ComputeErrorDomainFunctional(
    [[maybe_unused]] const GlobalVector& u,
    [[maybe_unused]] const DomainFunctional& F) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputeDomainFunctional\" not written!"
      << std::endl;
    abort();
  }

  virtual double ComputePointFunctional(
    [[maybe_unused]] const GlobalVector& u,
    [[maybe_unused]] const PointFunctional& FP) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputePointFunctional\" not written!"
      << std::endl;
    abort();
  }

  virtual void EvaluateCellRightHandSide(
    [[maybe_unused]] GlobalVector& f,
    [[maybe_unused]] const DomainRightHandSide& CF,
    [[maybe_unused]] double d = 1.) const
  {
    std::cerr
      << "\"DiscretizationInterface::EvaluateCellRighthandside\" not written!"
      << std::endl;
    abort();
  }

  virtual void EvaluateBoundaryCellRightHandSide(
    [[maybe_unused]] GlobalVector& f,
    [[maybe_unused]] const IntSet& Colors,
    [[maybe_unused]] const BoundaryRightHandSide& CF,
    [[maybe_unused]] double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::"
                 "EvaluateBoundaryCellRighthandside\" not written!"
              << std::endl;
    abort();
  }

  virtual void EvaluateParameterRightHandSide(
    [[maybe_unused]] GlobalVector& f,
    [[maybe_unused]] const DomainRightHandSide& CF,
    [[maybe_unused]] double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::EvaluateParameterRighthandside\" "
                 "not written!"
              << std::endl;
    abort();
  }

  virtual void EvaluateBoundaryParameterRightHandSide(
    [[maybe_unused]] GlobalVector& f,
    [[maybe_unused]] const IntSet& Colors,
    [[maybe_unused]] const BoundaryRightHandSide& CF,
    [[maybe_unused]] double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::"
                 "EvaluateBoundaryParameterRighthandside\" not written!"
              << std::endl;
    abort();
  }

  virtual void InterpolateDomainFunction(
    [[maybe_unused]] GlobalVector& f,
    [[maybe_unused]] const DomainFunction& DF) const
  {
    std::cerr
      << "\"DiscretizationInterface::InterpolateDomainFunction\" not written!"
      << std::endl;
    abort();
  }

  virtual void InterpolateCellDomainFunction(
    [[maybe_unused]] GlobalVector& f,
    [[maybe_unused]] const DomainFunction& DF) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateCellDomainFunction\" "
                 "not written!"
              << std::endl;
    abort();
  }

  virtual void ConstructInterpolator(
    [[maybe_unused]] MgInterpolatorInterface* I,
    [[maybe_unused]] const MeshTransferInterface* MT)
  {
    std::cerr
      << "\"DiscretizationInterface::ConstructInterpolator\" not written!"
      << std::endl;
    abort();
  }

  virtual void GetVolumes([[maybe_unused]] DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetVolumes\" not written!"
              << std::endl;
    abort();
  }

  virtual void GetAreas([[maybe_unused]] DoubleVector& a,
                        [[maybe_unused]] const IntSet& Colors) const
  {
    std::cerr << "\"DiscretizationInterface::GetAreas\" not written!"
              << std::endl;
    abort();
  }

  virtual void GetMassDiag([[maybe_unused]] DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetMassDiag\" not written!"
              << std::endl;
    abort();
  }

  virtual void GetBoundaryMassDiag([[maybe_unused]] DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetBoundaryMassDiag\" not written!"
              << std::endl;
    abort();
  }

  virtual void RhsCurve([[maybe_unused]] GlobalVector& F,
                        [[maybe_unused]] const Curve& C,
                        [[maybe_unused]] int comp,
                        [[maybe_unused]] int N) const
  {
    std::cerr << "\"DiscretizationInterface::RhsCurve\" not written!"
              << std::endl;
    abort();
  }
};
} // namespace Gascoigne

#endif
