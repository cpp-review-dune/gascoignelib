/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the Gascoigne 3D
 *authors
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

#ifndef __StdSolver_h
#define __StdSolver_h

#include "domainfunction.h"
#include "gascoigne.h"
#include "gascoignemesh.h"
#include "gascoignevisualization.h"

#include "ghostagent.h"
#include "iluagent.h"

#include "hierarchicalmesh.h"
#include "multigridmeshinterface.h"
#include "pointfunctional.h"
#include "pressurefilter.h"
#include "residualfunctional.h"
#include "solverdata.h"
#include "sparsestructure.h"
#include "stopwatch.h"

#include "discretizationinterface.h"
#include "iluinterface.h"
#include "matrixinterface.h"
#include "problemdescriptorinterface.h"
#include "vectorinterface.h"

/*-----------------------------------------*/

namespace Gascoigne {
//////////////////////////////////////////////
///
///@brief
/// Default nonlinear StdSolver

///
///
//////////////////////////////////////////////

class StdSolver
{
private:
  //
  //   Daten
  //

  // 0.

#ifdef __WITH_THREADS__
  IndexType __n_threads, __min_patches_per_thread;
  bool __with_thread_ilu;
  std::vector<std::vector<IndexType>> __thread_domain2node;
  // First component of pair is domain in wich the node lies, second is
  // local index of this node in __thread_domain2node
  std::vector<std::vector<std::pair<IndexType, IndexType>>>
    __thread_node2domain;
#endif

  // 1. Gitter

  const GascoigneMesh* _MP;
  const HierarchicalMesh* _HM;
  std::map<IndexType, IndexType> _PeriodicPairs;

#ifdef USE_CUDA
  static constexpr bool use_cuda = true;
#else
  static constexpr bool use_cuda = false;
#endif
  // 2. Matrizen

  //  MatrixInterface* _MAP;
  //  IluInterface* _MIP;

protected:
#ifdef __WITH_THREADS__
  int NThreads() const
  {
    return __n_threads;
  }
#endif

  // 3. Discretization

  DiscretizationInterface* _ZP;

  // 4. Vectors and Matrices
  mutable GhostVectorAgent vector_agent;
  mutable MatrixAgent matrix_agent;
  mutable IluAgent ilu_agent;

  // 5. Anwendungsklassen

  const ProblemDescriptorInterface* _PDX;
  //    const NumericInterface *_NI;

  // 6. Steuerparameter

  bool _distribute;

  mutable IndexType _ndirect;
  mutable bool _directsolver;
  mutable std::string _discname;
  mutable std::string _matrixtype;

  mutable IndexType _PrimalSolve;
  ParamFile _paramfile;

  bool _useUMFPACK;

  // 5. sonstiges

  PressureFilter _PF;
  /*   double               omega_domain; */

  //
  //        Funktionen
  //

  // 0. Zugriff

  const GascoigneMesh*& GetMeshPointer()
  {
    return _MP;
  }

  // virtual SolverData& GetSolverData()
  // {
  //   return _Dat;
  // }
  virtual const SolverData& GetSolverData() const
  {
    return GetProblemDescriptor()->GetSolverData();
  }
  virtual PressureFilter& GetPfilter()
  {
    return _PF;
  }
  virtual const PressureFilter& GetPfilter() const
  {
    return _PF;
  }

  // 0.3 Matrizen

  // MatrixInterface*& GetMatrixPointer()
  // {
  //   return _MAP;
  // }

  virtual DiscretizationInterface*& GetDiscretizationPointer()
  {
    return _ZP;
  }

  // virtual IluInterface*& GetIluPointer()
  // {
  //   return _MIP;
  // }

  // 1. Initialisierung

  virtual void SetDefaultValues(std::string discname,
                                std::string matrixtype,
                                IndexType ndirect);

  virtual DiscretizationInterface* NewDiscretization(
    IndexType dimension,
    const std::string& discname);
  virtual MatrixInterface* NewMatrix(IndexType ncomp,
                                     const std::string& matrixtype);
  virtual IluInterface* NewIlu(const Matrix& A,
                               IndexType ncomp,
                               const std::string& matrixtype);

  //
  /// new interface-function for individual size of vectors
  //

  virtual void smooth(IndexType niter,
                      const Matrix& A,
                      Vector& x,
                      const Vector& y,
                      Vector& h) const;
  virtual void PermutateIlu(Matrix& A, const Vector& gu) const;
  virtual void modify_ilu(IluInterface& I, IndexType ncomp) const;

  virtual DoubleVector IntegrateSolutionVector(const Vector& u) const;
  virtual void _check_consistency(const Equation* EQ,
                                  const DiscretizationInterface* DI) const;
  virtual void DirichletMatrixOnlyRow(Matrix& A) const;

public:
  StdSolver();
  virtual ~StdSolver();

  virtual std::string GetName() const
  {
    return "StdSolver";
  }

  virtual void BasicInit(const ParamFile& paramfile, const IndexType dimension);
  ////                const NumericInterface *NI);
  virtual void SetProblem(const ProblemDescriptorInterface& PDX);
  virtual void SetDiscretization(DiscretizationInterface& DI,
                                 bool init = false);
  virtual const ProblemDescriptorInterface* GetProblemDescriptor() const
  {
    assert(_PDX);
    if (_PDX) {
      return _PDX;
    }
    return nullptr;
  }
  virtual const ParamFile& GetParamfile() const
  {
    return _paramfile;
  }

  virtual void NewMesh(const GascoigneMesh* MP);

  virtual const GascoigneMesh* GetMesh() const
  {
    return _MP;
  }

  // 0.2 Discretization

  virtual const DiscretizationInterface* GetDiscretization() const
  {
    assert(_ZP);
    return _ZP;
  }
  virtual DiscretizationInterface* GetDiscretization()
  {
    assert(_ZP);
    return _ZP;
  }

  virtual void AddPeriodicNodes(SparseStructure* SA);

  // virtual IluInterface* GetIlu() const
  // {
  //   assert(_MIP);
  //   return _MIP;
  // }

  virtual bool DirectSolver() const
  {
    return _directsolver;
  }

  virtual void AddNodeVector(const std::string& name, const Vector& q)
  {
    assert(q.GetType() == "node");
    GetDiscretization()->AddNodeVector(name, &GetGV(q));
  }
  virtual void AddCellVector(const std::string& name, const Vector& q)
  {
    assert(q.GetType() == "cell");
    GetDiscretization()->AddCellVector(name, &GetGV(q));
  }
  virtual void AddParameterVector(const std::string& name,
                                  const GlobalParameterVector* q)
  {
    GetDiscretization()->AddParameterVector(name, q);
  }
  virtual void DeleteNodeVector(const std::string& name)
  {
    GetDiscretization()->DeleteNodeVector(name);
  }
  virtual void DeleteCellVector(const std::string& name)
  {
    GetDiscretization()->DeleteCellVector(name);
  }
  virtual void DeleteParameterVector(const std::string& name)
  {
    GetDiscretization()->DeleteParameterVector(name);
  }

  virtual void OutputSettings() const;
  virtual void PointVisu(const std::string& name,
                         const GlobalVector& u,
                         IndexType i) const;
  virtual void CellVisu(const std::string& name,
                        const GlobalVector& u,
                        IndexType i) const;

  virtual void ConstructInterpolator(MgInterpolatorInterface* I,
                                     const MeshTransferInterface* MT);
  virtual void VisuGrid(const std::string& name, IndexType i) const;

  //
  /// vector & Matrix  - manamgement
  //

  virtual void ReInitMatrix(const Matrix& A);

  virtual void RegisterVector(const Vector& g);
  virtual void ReInitVector(Vector& dst);
  virtual void ReInitVector(Vector& dst, IndexType comp);

  // Access to Vector & Matrix Data
  virtual GlobalVector& GetGV(Vector& u) const
  {
    return vector_agent(u);
  }
  virtual const GlobalVector& GetGV(const Vector& u) const
  {
    return vector_agent(u);
  }

  // Access to Vector & Matrix Data
  virtual MatrixInterface& GetMatrix(Matrix& A) const
  {
    return matrix_agent(A);
  }
  virtual const MatrixInterface& GetMatrix(const Matrix& A) const
  {
    return matrix_agent(A);
  }
  virtual IluInterface& GetIlu(Matrix& A) const
  {
    return ilu_agent(A);
  }
  virtual const IluInterface& GetIlu(const Matrix& A) const
  {
    return ilu_agent(A);
  }

  //
  /// vector - hanging nodes
  //

  virtual bool GetDistribute() const
  {
    return _distribute;
  }
  virtual void SetDistribute(bool dist)
  {
    _distribute = dist;
  }

  virtual void HNAverage(const Vector& x) const;
  virtual void HNZero(const Vector& x) const;
  virtual void HNDistribute(Vector& x) const;
  virtual void HNAverageData() const;
  virtual void HNZeroData() const;

  //
  /// vector - io
  //

  virtual void Visu(const std::string& name,
                    const Vector& u,
                    IndexType i) const;
  virtual void Write(const Vector& u, const std::string& filename) const;
  virtual void Read(Vector& u, const std::string& filename) const;

  //
  /// vector - interpolation
  //

  virtual void InterpolateSolution(Vector& u, const GlobalVector& uold) const;

  //
  /// vector - rhs (integration)
  //

  virtual void Rhs(Vector& f, double d = 1.) const;

  //
  /// vector - residual (integration)
  //

  virtual void Form(Vector& y, const Vector& x, double d) const;
  virtual void AdjointForm(Vector& y, const Vector& x, double d) const;

  //
  /// vector - boundary condition
  //
  virtual void SetBoundaryVector(Vector& f) const;
  virtual void SetPeriodicVector(Vector& f) const;
  virtual void SetBoundaryVectorZero(Vector& Gf) const;
  virtual void SetBoundaryVectorStrong(Vector& f,
                                       const BoundaryManager& BM,
                                       const DirichletData& DD,
                                       double d = 1.) const;
  virtual void SetPeriodicVectorStrong(Vector& f,
                                       const BoundaryManager& BM,
                                       const PeriodicData& PD,
                                       double d = 1.) const;
  virtual void SetPeriodicVectorZero(Vector& gf) const;

  //
  /// vector - linear algebra
  //

  virtual double NewtonNorm(const Vector& u) const;
  virtual void residualgmres(const Matrix& A,
                             Vector& y,
                             const Vector& x,
                             const Vector& b) const;
  virtual void MatrixResidual(const Matrix& A,
                              Vector& y,
                              const Vector& x,
                              const Vector& b) const;
  virtual void vmult(const Matrix& A,
                     Vector& y,
                     const Vector& x,
                     double d) const;
  virtual void vmulteq(const Matrix& A,
                       Vector& y,
                       const Vector& x,
                       double d) const;
  virtual void smooth_pre(const Matrix& A,
                          Vector& y,
                          const Vector& x,
                          Vector& h) const;
  virtual void smooth_exact(const Matrix& A,
                            Vector& y,
                            const Vector& x,
                            Vector& h) const;
  virtual void smooth_post(const Matrix& A,
                           Vector& y,
                           const Vector& x,
                           Vector& h) const;
  virtual void Zero(Vector& dst) const;

  //
  /// vector - additional
  //

  virtual void SubtractMean(Vector& x) const;
  virtual void SubtractMeanAlgebraic(Vector& x) const;

  //
  /// vector - matrix
  //

  virtual void AssembleMatrix(Matrix& A, const Vector& u, double d) const;
  virtual void DirichletMatrix(Matrix& A) const;
  virtual void PeriodicMatrix(Matrix& A) const;
  virtual void MatrixZero(Matrix& A) const;
  virtual void ComputeIlu(Matrix& A, const Vector& u) const;
  virtual void AssembleDualMatrix(Matrix& A, const Vector& gu, double d);
  virtual void MassMatrixVector(Vector& f, const Vector& gu, double d) const
  {
    abort();
  }

  // virtual MatrixInterface* GetMatrix() const
  // {
  //   return _MAP;
  // }

  //
  /// vector - "postprocessing"
  //

  virtual void ComputeError(const Vector& u, GlobalVector& err) const;
  virtual void AssembleError(GlobalVector& eta,
                             const Vector& u,
                             GlobalVector& err) const;
  virtual double ComputeFunctional(Vector& f,
                                   const Vector& u,
                                   const Functional* FP);

  virtual double ComputeBoundaryFunctional(Vector& f,
                                           const Vector& u,
                                           Vector& z,
                                           const BoundaryFunctional* FP) const;
  virtual double ComputeDomainFunctional(const Vector& u,
                                         const DomainFunctional* FP) const;
  virtual double ComputePointFunctional(Vector& f,
                                        const Vector& u,
                                        Vector& z,
                                        const PointFunctional* NFP) const;
  virtual double ComputeResidualFunctional(Vector& f,
                                           const Vector& u,
                                           Vector& z,
                                           const ResidualFunctional* FP) const;
  virtual void EvaluateCellRightHandSide(Vector& f,
                                         const DomainRightHandSide& CF,
                                         double d = 1.) const;
  virtual void EvaluateBoundaryCellRightHandSide(
    Vector& f,
    const BoundaryRightHandSide& CF,
    const BoundaryManager& BM,
    double d = 1.) const;
  virtual void EvaluateParameterRightHandSide(Vector& f,
                                              const DomainRightHandSide& CF,
                                              double d = 1.) const;
  virtual void EvaluateBoundaryParameterRightHandSide(
    Vector& f,
    const BoundaryRightHandSide& CF,
    const BoundaryManager& BM,
    double d = 1.) const;
  virtual void InterpolateDomainFunction(Vector& f,
                                         const DomainFunction& DF) const;

  //
  /// vector - initialize
  //

  virtual void BoundaryInit(Vector& u) const;
  virtual void SolutionInit(Vector& u) const;

  //
  /// HierarchicalMesh
  //

  virtual const HierarchicalMesh*& GetHierarchicalMeshPointer()
  {
    return _HM;
  }
  virtual const HierarchicalMesh* GetHierarchicalMesh() const
  {
    return _HM;
  }

  //
  /// for gmres
  //
  virtual void DeleteVector(Vector& p) const;

  virtual double ScalarProduct(const Vector& y, const Vector& x) const;
  virtual void Equ(Vector& dst, double s, const Vector& src) const;
  virtual void Add(Vector& dst, double s, const Vector& src) const;
  virtual void SAdd(double s1, Vector& dst, double s2, const Vector& src) const;
  virtual double Norm(const Vector& dst) const;

  virtual void RhsCurve(Vector& f,
                        const Curve& C,
                        IndexType comp,
                        IndexType N) const;
  virtual double ScalarProductWithFluctuations(DoubleVector& eta,
                                               const Vector& gf,
                                               const Vector& gz) const;

#ifdef __WITH_THREADS__
  //
  /// For Threads
  //
  virtual void ThreadPartitionMesh();
#endif
};
} // namespace Gascoigne

#endif
