/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the Gascoigne
 *3D authors
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

#ifndef __SolverInterface_h
#define __SolverInterface_h

#include "gascoigne.h"
#include "iluinterface.h"
#include "matrixinterface.h"
#include "sparsestructure.h"

#include "discretizationinterface.h"
#include "functional.h"
#include "multigridmeshinterface.h"

#include "curve.h"
#include "meshinterface.h"
#include "meshtransferinterface.h"
#include "mginterpolatorinterface.h"
//#include "numericinterface.h"
#include "paramfile.h"
#include "problemdescriptorinterface.h"
#include "vectorinterface.h"

/*---------------------------------------------------------*/

namespace Gascoigne {

//////////////////////////////////////////////
///
///@brief
/// Interface class for Solver

///  Some porporties
///  - lives on one level of the hierarchy
///  - stores the matrices
///  - holds memory for vectors
///  - provides e.g. nonlinear and linear residuals of the equations
///  - calls class Discretization
///
//////////////////////////////////////////////

class SolverInterface
{
private:
protected:
public:
  SolverInterface() {}
  virtual ~SolverInterface() {}

  virtual std::string GetName() const = 0;

  virtual void BasicInit(const ParamFile* paramfile,
                         const int dimension,
                         const NumericInterface* NI = NULL) = 0;

  virtual void SetProblem(const ProblemDescriptorInterface& PD) = 0;
  virtual void SetDiscretization(DiscretizationInterface& DI,
                                 bool init = false) = 0;
  virtual const ProblemDescriptorInterface* GetProblemDescriptor() const = 0;
  virtual const ParamFile* GetParamfile() const = 0;

  virtual bool DirectSolver() const = 0;

  virtual void NewMesh(const MeshInterface* MP) = 0;
  virtual const MeshInterface* GetMesh() const = 0;

  virtual void RegisterMatrix() = 0;
  virtual void ReInitMatrix() = 0;

  virtual void OutputSettings() const = 0;
  virtual void ConstructInterpolator(MgInterpolatorInterface* I,
                                     const MeshTransferInterface* MT) = 0;

  virtual void VisuGrid(const std::string& name, int i) const = 0;
  virtual void SetTimeData(double dt,
                           double theta,
                           double time,
                           double oldrhs = -1,
                           double newrhs = -1){};

  virtual void AddNodeVector(const std::string&, const Vector& q) = 0;
  virtual void AddCellVector(const std::string&, const Vector& q) = 0;
  virtual void AddParameterVector(const std::string&,
                                  const GlobalParameterVector* q) = 0;
  virtual void DeleteNodeVector(const std::string&) = 0;
  virtual void DeleteCellVector(const std::string&) = 0;
  virtual void DeleteParameterVector(const std::string&) = 0;
  virtual void DeleteVector(Gascoigne::Vector& p) const = 0;

  virtual const DiscretizationInterface* GetDiscretization() const = 0;
  virtual DiscretizationInterface* GetDiscretization() = 0;

  //
  /// vector - manamgement
  //
  virtual void RegisterVector(const Vector& g) = 0;
  virtual GlobalVector& GetGV(Vector& u) const = 0;
  virtual const GlobalVector& GetGV(const Vector& u) const = 0;
  virtual void ReInitVector(Vector& dst) = 0;
  virtual void ReInitVector(Vector& dst, int comp) = 0;

  //
  /// vector - hanging nodes
  //
  virtual bool GetDistribute() const { abort(); }
  virtual void SetDistribute(bool dist) { abort(); }

  virtual void HNAverage(const Vector& x) const = 0;
  virtual void HNZero(const Vector& x) const = 0;
  virtual void HNDistribute(Vector& x) const = 0;

  virtual void HNAverageData() const = 0;
  virtual void HNZeroData() const = 0;

  //
  /// vector - io
  //
  virtual void Visu(const std::string& name, const Vector& u, int i) const = 0;
  virtual void Write(const Vector& u, const std::string& filename) const = 0;
  virtual void Read(Vector& u, const std::string& filename) const = 0;

  //
  /// vector - interpolation
  //
  virtual void InterpolateSolution(Vector& u,
                                   const GlobalVector& uold) const = 0;

  //
  /// vector - rhs (integration)
  //
  virtual void Rhs(Vector& f, double d = 1.) const = 0;
  virtual void TimeRhsOperator(Vector& f, const Vector& u) const
  {
    std::cerr << "\"SolverInterface::TimeRhsOperator\" not written!"
              << std::endl;
    abort();
  }
  virtual void TimeRhs(int k, Vector& f) const
  {
    std::cerr << "\"SolverInterface::TimeRhs\" not written!" << std::endl;
    abort();
  }
  virtual void InitialCondition(Vector& f, double d = 1.) const
  {
    std::cerr << "\"SolverInterface::IC\" not written!" << std::endl;
    abort();
  }
  virtual void L2Projection(Vector& u, Vector& f)
  {
    std::cerr << "\"SolverInterface::L2Projection\" not written!" << std::endl;
    abort();
  }

  virtual void RhsCurve(Vector& f, const Curve& C, int comp, int N) const
  {
    std::cerr << "\"SolverInterface::RhsCurve\" not written!" << std::endl;
    abort();
  }

  //
  /// vector - residual (integration)
  //
  virtual void Form(Vector& y, const Vector& x, double d) const = 0;
  virtual void AdjointForm(Vector& y, const Vector& x, double d) const
  {
    std::cerr << "\"SolverInterface::AdjointForm\" not written!" << std::endl;
    abort();
  }

  //
  /// vector - boundary condition
  //
  virtual void SetBoundaryVector(Vector& f) const = 0;
  virtual void SetPeriodicVector(Vector& f) const = 0;
  virtual void SetBoundaryVectorZero(Vector& f) const = 0;
  virtual void SetBoundaryVectorStrong(Vector& f,
                                       const BoundaryManager& BM,
                                       const DirichletData& DD,
                                       double d = 1.) const = 0;
  virtual void SetPeriodicVectorStrong(Vector& f,
                                       const BoundaryManager& BM,
                                       const PeriodicData& PD,
                                       double d = 1.) const = 0;
  virtual void SetPeriodicVectorZero(Vector& f) const = 0;
  virtual void AddPeriodicNodes(SparseStructure* SA)
  {
    std::cerr << "\"SolverInterface::AddPeriodicNodes\" not written!"
              << std::endl;
    abort();
  }

  //
  /// vector - linear algebra
  //
  virtual double NewtonNorm(const Vector& u) const = 0;

  virtual void residualgmres(Vector& y,
                             const Vector& x,
                             const Vector& b) const = 0;
  virtual void MatrixResidual(Vector& y,
                              const Vector& x,
                              const Vector& b) const = 0;
  virtual void vmult(Vector& y, const Vector& x, double d) const = 0;
  virtual void vmulteq(Vector& y, const Vector& x, double d) const = 0;
  virtual void smooth(int niter,
                      Vector& y,
                      const Vector& x,
                      Vector& h) const = 0;
  virtual void smooth_pre(Vector& y, const Vector& x, Vector& h) const = 0;
  virtual void smooth_exact(Vector& y, const Vector& x, Vector& h) const = 0;
  virtual void smooth_post(Vector& y, const Vector& x, Vector& h) const = 0;
  virtual void Zero(Vector& dst) const = 0;

  //
  /// vector - additional
  //
  virtual void SubtractMean(Vector& x) const = 0;
  virtual void SubtractMeanAlgebraic(Vector& x) const = 0;

  //
  /// vector - matrix
  //
  virtual void AssembleMatrix(const Vector& u, double d = 1.) = 0;
  virtual void DirichletMatrix() const = 0;
  virtual void MatrixZero() const = 0;
  virtual void ComputeIlu(const Vector& u) const = 0;
  virtual void ComputeIlu() const = 0;
  virtual void AssembleDualMatrix(const Vector& gu, double d) = 0;
  virtual void MassMatrixVector(Vector& f,
                                const Vector& gu,
                                double d) const = 0;
  virtual void InverseMassMatrix(Vector& u, const Vector& f) const {};
  virtual void PeriodicMatrix() const = 0;

  //
  /// vector
  //
  virtual double ScalarProduct(const Vector& y, const Vector& x) const = 0;
  virtual void Equ(Vector& dst, double s, const Vector& src) const = 0;
  virtual void Add(Vector& dst, double s, const Vector& src) const = 0;
  virtual void SAdd(double s1,
                    Vector& dst,
                    double s2,
                    const Vector& src) const = 0;
  virtual double Norm(const Vector& dst) const = 0;

  //
  /// vector - "postprocessing"
  //
  virtual void ComputeError(const Vector& u, GlobalVector& err) const = 0;
  virtual void AssembleError(GlobalVector& eta,
                             const Vector& u,
                             GlobalVector& err) const = 0;
  virtual double ComputeFunctional(Vector& f,
                                   const Vector& u,
                                   const Functional* FP) = 0;

  //
  /// vector - initialize
  //
  virtual void BoundaryInit(Vector& u) const = 0;
  virtual void SolutionInit(Vector& u) const = 0;

  virtual double ScalarProductWithFluctuations(DoubleVector& eta,
                                               const Vector& gf,
                                               const Vector& gz) const = 0;
};
} // namespace Gascoigne
#endif
