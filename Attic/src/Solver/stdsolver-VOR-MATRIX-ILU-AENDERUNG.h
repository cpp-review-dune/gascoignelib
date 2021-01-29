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
#include "ghostvectoragent.h"
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

/*-----------------------------------------*/

namespace Gascoigne
{
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
  int __n_threads, __min_patches_per_thread;
  bool __with_thread_ilu;
  std::vector<std::vector<int>> __thread_domain2node;
  // First component of pair is domain in wich the node lies, second is
  // local index of this node in __thread_domain2node
  std::vector<std::vector<std::pair<int, int>>> __thread_node2domain;
#endif

  // 1. Gitter

  const GascoigneMesh* _MP;
  const HierarchicalMesh* _HM;
  std::map<int, int> _PeriodicPairs;

#ifdef USE_CUDA
  static constexpr bool use_cuda = true;
#else
  static constexpr bool use_cuda = false;
#endif
  // 2. Matrizen

  MatrixInterface* _MAP;
  IluInterface* _MIP;

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

  mutable GhostVectorAgent _NGVA;
  mutable MatrixAgent      matrix_agent;

  // 5. Anwendungsklassen

  const ProblemDescriptorInterface* _PDX;
  //    const NumericInterface *_NI;

  // 6. Steuerparameter

  bool _distribute;

  mutable int _ndirect;
  mutable bool _directsolver;
  mutable std::string _discname;
  mutable std::string _matrixtype;

  SolverData _Dat;
  mutable int _PrimalSolve;
  const ParamFile* _paramfile;

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

  virtual SolverData& GetSolverData()
  {
    return _Dat;
  }
  virtual const SolverData& GetSolverData() const
  {
    return _Dat;
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

  MatrixInterface*& GetMatrixPointer()
  {
    return _MAP;
  }

  virtual DiscretizationInterface*& GetDiscretizationPointer()
  {
    return _ZP;
  }

  virtual IluInterface*& GetIluPointer()
  {
    return _MIP;
  }

  // 1. Initialisierung

  virtual void SetDefaultValues(std::string discname, std::string matrixtype,
                                int ndirect);

  virtual DiscretizationInterface* NewDiscretization(int dimension,
                                                     const std::string& discname);
  virtual MatrixInterface* NewMatrix(int ncomp, const std::string& matrixtype);
  virtual IluInterface* NewIlu(int ncomp, const std::string& matrixtype);

  virtual void RegisterMatrix(int ncomp);

  //
  /// new interface-function for individual size of vectors
  //

  virtual void smooth(int niter, VectorInterface& x, const VectorInterface& y,
                      VectorInterface& h) const;
  virtual void PermutateIlu(const VectorInterface& gu) const;
  virtual void modify_ilu(IluInterface& I, int ncomp) const;

  virtual DoubleVector IntegrateSolutionVector(const VectorInterface& u) const;
  virtual void _check_consistency(const Equation* EQ,
                                  const DiscretizationInterface* DI) const;
  virtual void DirichletMatrixOnlyRow() const;

public:
  StdSolver();
  virtual ~StdSolver();

  virtual std::string GetName() const
  {
    return "StdSolver";
  }

  virtual void BasicInit(const ParamFile* paramfile, const int dimension);
  ////                const NumericInterface *NI);
  virtual void SetProblem(const ProblemDescriptorInterface& PDX);
  virtual void SetDiscretization(DiscretizationInterface& DI, bool init = false);
  virtual const ProblemDescriptorInterface* GetProblemDescriptor() const
  {
    assert(_PDX);
    return _PDX;
  }
  virtual const ParamFile* GetParamfile() const
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

  virtual void ReInitMatrix();
  virtual void AddPeriodicNodes(SparseStructure* SA);

  virtual IluInterface* GetIlu() const
  {
    assert(_MIP);
    return _MIP;
  }

   virtual bool DirectSolver() const
  {
    return _directsolver;
  }

  virtual void AddNodeVector(const std::string& name, const VectorInterface& q)
  {
    assert(q.GetType() == "node");
    GetDiscretization()->AddNodeVector(name, &GetGV(q));
  }
  virtual void AddCellVector(const std::string& name, const VectorInterface& q)
  {
    assert(q.GetType() == "cell");
    GetDiscretization()->AddCellVector(name, &GetGV(q));
  }
  virtual void AddParameterVector(const std::string& name, const GlobalParameterVector* q)
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
  virtual void PointVisu(const std::string& name, const GlobalVector& u, int i) const;
  virtual void CellVisu(const std::string& name, const GlobalVector& u, int i) const;

  virtual void ConstructInterpolator(MgInterpolatorInterface* I,
                                     const MeshTransferInterface* MT);
  virtual void VisuGrid(const std::string& name, int i) const;

  //
  /// vector & Matrix  - manamgement
  //

  virtual void RegisterMatrix();

  virtual void ReInitMatrix(Matrix& A);
  
  virtual void RegisterVector(const VectorInterface& g);
  virtual void ReInitVector(VectorInterface& dst);
  virtual void ReInitVector(VectorInterface& dst, int comp);

  virtual GlobalVector& GetGV(VectorInterface& u) const
  {
    return _NGVA(u);
  }
  virtual const GlobalVector& GetGV(const VectorInterface& u) const
  {
    return _NGVA(u);
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

  virtual void HNAverage(const VectorInterface& x) const;
  virtual void HNZero(const VectorInterface& x) const;
  virtual void HNDistribute(VectorInterface& x) const;
  virtual void HNAverageData() const;
  virtual void HNZeroData() const;

  //
  /// vector - io
  //

  virtual void Visu(const std::string& name, const VectorInterface& u, int i) const;
  virtual void Write(const VectorInterface& u, const std::string& filename) const;
  virtual void Read(VectorInterface& u, const std::string& filename) const;

  //
  /// vector - interpolation
  //

  virtual void InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const;

  //
  /// vector - rhs (integration)
  //

  virtual void Rhs(VectorInterface& f, double d = 1.) const;

  //
  /// vector - residual (integration)
  //

  virtual void Form(VectorInterface& y, const VectorInterface& x, double d) const;
  virtual void AdjointForm(VectorInterface& y, const VectorInterface& x, double d) const;

  //
  /// vector - boundary condition
  //
  virtual void SetBoundaryVector(VectorInterface& f) const;
  virtual void SetPeriodicVector(VectorInterface& f) const;
  virtual void SetBoundaryVectorZero(VectorInterface& Gf) const;
  virtual void SetBoundaryVectorStrong(VectorInterface& f, const BoundaryManager& BM,
                                       const DirichletData& DD, double d = 1.) const;
  virtual void SetPeriodicVectorStrong(VectorInterface& f, const BoundaryManager& BM,
                                       const PeriodicData& PD, double d = 1.) const;
  virtual void SetPeriodicVectorZero(VectorInterface& gf) const;

  //
  /// vector - linear algebra
  //

  virtual double NewtonNorm(const VectorInterface& u) const;
  virtual void residualgmres(VectorInterface& y, const VectorInterface& x,
                             const VectorInterface& b) const;
  virtual void MatrixResidual(VectorInterface& y, const VectorInterface& x,
                              const VectorInterface& b) const;
  virtual void vmult(VectorInterface& y, const VectorInterface& x, double d) const;
  virtual void vmulteq(VectorInterface& y, const VectorInterface& x, double d) const;
  virtual void smooth_pre(VectorInterface& y, const VectorInterface& x,
                          VectorInterface& h) const;
  virtual void smooth_exact(VectorInterface& y, const VectorInterface& x,
                            VectorInterface& h) const;
  virtual void smooth_post(VectorInterface& y, const VectorInterface& x,
                           VectorInterface& h) const;
  virtual void Zero(VectorInterface& dst) const;

  //
  /// vector - additional
  //

  virtual void SubtractMean(VectorInterface& x) const;
  virtual void SubtractMeanAlgebraic(VectorInterface& x) const;

  //
  /// vector - matrix
  //

  virtual void AssembleMatrix(const VectorInterface& u, double d);
  virtual void DirichletMatrix() const;
  virtual void PeriodicMatrix() const;
  virtual void MatrixZero() const;
  virtual void ComputeIlu(const VectorInterface& u) const;
  virtual void ComputeIlu() const;
  virtual void AssembleDualMatrix(const VectorInterface& gu, double d);
  virtual void MassMatrixVector(VectorInterface& f, const VectorInterface& gu,
                                double d) const
  {
    abort();
  }

  virtual MatrixInterface* GetMatrix() const
  {
    return _MAP;
  }

  //
  /// vector - "postprocessing"
  //

  virtual void ComputeError(const VectorInterface& u, GlobalVector& err) const;
  virtual void AssembleError(GlobalVector& eta, const VectorInterface& u,
                             GlobalVector& err) const;
  virtual double ComputeFunctional(VectorInterface& f, const VectorInterface& u,
                                   const Functional* FP);

  virtual double ComputeBoundaryFunctional(VectorInterface& f, const VectorInterface& u,
                                           VectorInterface& z,
                                           const BoundaryFunctional* FP) const;
  virtual double ComputeDomainFunctional(const VectorInterface& u,
                                         const DomainFunctional* FP) const;
  virtual double ComputePointFunctional(VectorInterface& f, const VectorInterface& u,
                                        VectorInterface& z,
                                        const PointFunctional* NFP) const;
  virtual double ComputeResidualFunctional(VectorInterface& f, const VectorInterface& u,
                                           VectorInterface& z,
                                           const ResidualFunctional* FP) const;
  virtual void EvaluateCellRightHandSide(VectorInterface& f,
                                         const DomainRightHandSide& CF,
                                         double d = 1.) const;
  virtual void EvaluateBoundaryCellRightHandSide(VectorInterface& f,
                                                 const BoundaryRightHandSide& CF,
                                                 const BoundaryManager& BM,
                                                 double d = 1.) const;
  virtual void EvaluateParameterRightHandSide(VectorInterface& f,
                                              const DomainRightHandSide& CF,
                                              double d = 1.) const;
  virtual void EvaluateBoundaryParameterRightHandSide(VectorInterface& f,
                                                      const BoundaryRightHandSide& CF,
                                                      const BoundaryManager& BM,
                                                      double d = 1.) const;
  virtual void InterpolateDomainFunction(VectorInterface& f,
                                         const DomainFunction& DF) const;

  //
  /// vector - initialize
  //

  virtual void BoundaryInit(VectorInterface& u) const;
  virtual void SolutionInit(VectorInterface& u) const;

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
  virtual void DeleteVector(VectorInterface& p) const;

  virtual double ScalarProduct(const VectorInterface& y, const VectorInterface& x) const;
  virtual void Equ(VectorInterface& dst, double s, const VectorInterface& src) const;
  virtual void Add(VectorInterface& dst, double s, const VectorInterface& src) const;
  virtual void SAdd(double s1, VectorInterface& dst, double s2,
                    const VectorInterface& src) const;
  virtual double Norm(const VectorInterface& dst) const;

  virtual void RhsCurve(VectorInterface& f, const Curve& C, int comp, int N) const;
  virtual double ScalarProductWithFluctuations(DoubleVector& eta,
                                               const VectorInterface& gf,
                                               const VectorInterface& gz) const;

#ifdef __WITH_THREADS__
  //
  /// For Threads
  //
  virtual void ThreadPartitionMesh();
#endif
};
}  // namespace Gascoigne

#endif
