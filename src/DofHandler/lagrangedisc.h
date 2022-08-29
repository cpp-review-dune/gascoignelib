/*----------------------------   lagrangedisc.h ---------------------------*/
/*      $Id:$                 */
#ifndef __lagrangedisc_H
#define __lagrangedisc_H
/*----------------------------   lagrangedisc.h ---------------------------*/

/**
 *
 * Copyright (C) 2019,2020 by the Gascoigne 3D authors
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

#include "discretizationinterface.h"
#include "mginterpolatornested.h"
//#include "omp.h"
#include "cgbase.h"
#include "cgdofhandler.h"
#include "gascoignevisualization.h"
#include "pressurefilter.h"
#include "problemdescriptorbase.h"
#include "sparsestructure.h"
#include "stopwatch.h"

namespace Gascoigne {

/////////////////////////////////////////////
///
///@brief
///  ... comments DiscretizationInterface

///
///
/////////////////////////////////////////////

/**
 * General discretization based on the CGDofHandler
 * Continuous finite elements of dimension DIM
 * Each element has M^DIM dof's, e.g. local degree M-1
 **/
template<int DIM, int M, class FINITEELEMENT, class INTEGRATOR>
class LagrangeDisc : public DiscretizationInterface
{
private:
  CGDofHandler<DIM, M> dofhandler;
  mutable DataContainer datacontainer;

protected:
public:
  LagrangeDisc() {}
  ~LagrangeDisc() {}

  const CGDofHandler<DIM, M>& GetDofHandler() const { return dofhandler; }

  const DataContainer& GetDataContainer() const { return datacontainer; }
  void SetDataContainer(const DataContainer& q) { datacontainer = q; }

  //
  //// Functions called from the Solver
  //
  std::string GetName() const { return "Lagrange Discretization"; }

  void AddNodeVector(const std::string& name, const GlobalVector* q) const
  {
    datacontainer.AddNodeVector(name, q);
  }
  void DeleteNodeVector(const std::string& name) const
  {
    datacontainer.DeleteNodeVector(name);
  }

  void AddCellVector(const std::string& name, const GlobalVector* q) const
  {
    datacontainer.AddCellVector(name, q);
  }
  void DeleteCellVector(const std::string& name) const
  {
    datacontainer.DeleteCellVector(name);
  }

  void AddParameterVector(const std::string& name,
                          const GlobalParameterVector* q) const
  {
    datacontainer.AddParameterVector(name, q);
  }
  void DeleteParameterVector(const std::string& name) const
  {
    datacontainer.DeleteParameterVector(name);
  }

  void BasicInit(const ParamFile& pf)
  {
    // HANGING NODES
  }
  void ReInit(const GascoigneMesh* GM)
  {
    const DofHandler<DIM>* DH = dynamic_cast<const DofHandler<DIM>*>(GM);
    assert(DH);
    dofhandler.InitFromGascoigneMesh(*DH);
  }

  IndexType ndofs() const { return GetDofHandler().ndofs(); }
  IndexType nelements() const { return GetDofHandler().nelements(); }
  IndexType ndofs_withouthanging() const
  {
    abort();
    return ndofs();
    // HANGING NODES
  }
  Vertex2d vertex2d(IndexType i) const { return GetDofHandler().vertex2d(i); }
  Vertex3d vertex3d(IndexType i) const { return GetDofHandler().vertex3d(i); }

  // Visualization
  void VisuVtk(const ComponentInformation* CI,
               const ParamFile& pf,
               const std::string& name,
               const GlobalVector& u,
               IndexType i) const
  {
    std::cerr << "LagrangeDisc VisuVtk missing" << std::endl;

    // HNAverage(const_cast<GlobalVector&>(u));

    // GascoigneVisualization Visu;
    // Visu.SetMesh(GetDofHandler());
    // if (CI)
    // 	{
    // 	  Visu.AddPointVector(CI, &u);
    // 	}
    // else
    // 	{
    // 	  Visu.AddPointVector(&u);
    // 	}

    // Visu.read_parameters(&pf);
    // Visu.set_name(name);
    // Visu.step(i);
    // Visu.write();

    // HNZero(const_cast<GlobalVector&>(u));
  }

  //////////////////////////////////////////////////
  void Transformation(nmatrix<double>& T, IndexType iq) const
  {
    assert(GetDofHandler().dimension() == DIM);
    IndexType ne = GetDofHandler().dofs_per_element();

    IndexVector indices = GetDofHandler().GetElement(iq);
    assert(ne == indices.size());

    T.memory(DIM, ne);
    if (DIM == 2) {
      for (IndexType ii = 0; ii < ne; ii++) {
        Vertex2d v = GetDofHandler().vertex2d(indices[ii]);
        T(0, ii) = v.x();
        T(1, ii) = v.y();
      }
    } else if (DIM == 3) {
      for (IndexType ii = 0; ii < ne; ii++) {
        Vertex3d v = GetDofHandler().vertex3d(indices[ii]);
        T(0, ii) = v.x();
        T(1, ii) = v.y();
        T(2, ii) = v.z();
      }
    }
  }
  void ConstructInterpolator(MgInterpolatorInterface* I,
                             const MeshTransferInterface* MT)
  {
    return;
    abort();
    MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
    assert(IP);
    IP->BasicInit(MT);
  }

  void Structure(SparseStructureInterface* SI) const
  {
    SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
    assert(S);

    S->build_begin(ndofs());
    for (IndexType iq = 0; iq < GetDofHandler().nelements(); iq++) {
      IndexVector indices = GetDofHandler().GetElement(iq);
      // HANGING NODES
      // HN->CondenseHanging(indices);
      S->build_add(indices.begin(), indices.end());
    }
    S->build_end();
    // HANGING NODES
    // HN->SparseStructureDiag(S);
  }

  ////////////////////////////////////////////////// handling local / global
  void GlobalToGlobalData(LocalParameterData& QP) const
  {
    const GlobalParameterData& gpd = GetDataContainer().GetParameterData();
    QP.clear();

    for (auto p : gpd)
      QP.insert(make_pair(p.first, *p.second));
  }
  void LocalToGlobal(MatrixInterface& A,
                     EntryMatrix& E,
                     IndexType iq,
                     double s) const
  {
    IndexVector indices = GetDofHandler().GetElement(iq);

    // HANGING NODES
    // HN->CondenseHanging(E,indices);
    IndexVector::const_iterator start = indices.begin();
    IndexVector::const_iterator stop = indices.end();
#pragma omp critical
    A.entry(start, stop, E, s);
  }
  void LocalToGlobal(GlobalVector& f,
                     const LocalVector& F,
                     IndexType iq,
                     double s) const
  {
    IndexVector indices = GetDofHandler().GetElement(iq);
    for (IndexType ii = 0; ii < indices.size(); ii++) {
      IndexType i = indices[ii];
#pragma omp critical
      f.add_node(i, s, ii, F);
    }
  }
  void GlobalToLocal(LocalVector& U, const GlobalVector& u, IndexType iq) const
  {
    IndexVector indices = GetDofHandler().GetElement(iq);
    U.ReInit(u.ncomp(), indices.size());
    for (IndexType ii = 0; ii < indices.size(); ii++) {
      IndexType i = indices[ii];
      U.equ_node(ii, i, u);
    }
  }
  void GlobalToLocalCell(LocalVector& U,
                         const GlobalVector& u,
                         IndexType iq) const
  {
    U.ReInit(u.ncomp(), 1);
    for (IndexType c = 0; c < u.ncomp(); ++c)
      U(0, c) = u(iq, c);
  }
  void GlobalToLocalData(IndexType iq, LocalData& QN, LocalData& QC) const
  {
    const GlobalData& gnd = GetDataContainer().GetNodeData();
    QN.clear();
    GlobalData::const_iterator p = gnd.begin();
    for (; p != gnd.end(); p++) {
      GlobalToLocal(QN[p->first], *p->second, iq);
    }

    const GlobalData& gcd = GetDataContainer().GetCellData();
    QC.clear();
    GlobalData::const_iterator q = gcd.begin();
    for (; q != gcd.end(); q++) {
      GlobalToLocalCell(QC[q->first], *q->second, iq);
    }
  }
  void InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
  {
    abort();
  }

  ////////////////////////////////////////////////// integration ueber Zellen
  void Form(GlobalVector& f,
            const GlobalVector& u,
            const Equation& EQ,
            double d) const
  {
    LocalParameterData QP;
    GlobalToGlobalData(QP);

#pragma omp parallel
    {
      nmatrix<double> T;
      FINITEELEMENT finiteelement;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __U, __F;
      LocalData __QN, __QC;
      const auto EQCP = EQ.createNew();
      EQCP->SetParameterData(QP);

#pragma omp for schedule(static)
      for (IndexType iq = 0; iq < GetDofHandler().nelements(); ++iq) {
        Transformation(T, iq);

        finiteelement.ReInit(T);
        GlobalToLocal(__U, u, iq);
        GlobalToLocalData(iq, __QN, __QC);

        EQCP->point_cell(GetDofHandler().material(iq));

        integrator.Form(*EQCP, __F, finiteelement, __U, __QN, __QC);
        LocalToGlobal(f, __F, iq, d);
      }
      delete EQCP;
    }
  }
  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
  {
    LocalParameterData QP;
    GlobalToGlobalData(QP);

#pragma omp parallel
    {
      nmatrix<double> T;
      FINITEELEMENT finiteelement;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __F;
      LocalData __QN, __QC;

      auto* RHSCP = RHS.createNew();
      RHSCP->SetParameterData(QP);

#pragma omp for schedule(static)
      for (IndexType iq = 0; iq < GetDofHandler().nelements(); ++iq) {
        Transformation(T, iq);
        finiteelement.ReInit(T);

        GlobalToLocalData(iq, __QN, __QC);
        RHSCP->point_cell(GetDofHandler().material(iq));
        integrator.Rhs(*RHSCP, __F, finiteelement, __QN, __QC);
        LocalToGlobal(f, __F, iq, s);
      }
      delete RHSCP;
    }
  }

  void Matrix(MatrixInterface& A,
              const GlobalVector& u,
              const Equation& EQ,
              double d) const
  {
    LocalParameterData QP;
    GlobalToGlobalData(QP);

#pragma omp parallel
    {
      nmatrix<double> T;
      FINITEELEMENT finiteelement;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __U;
      LocalData __QN, __QC;
      EntryMatrix __E;

      auto* EQCP = EQ.createNew();
      EQCP->SetParameterData(QP);

#pragma omp for schedule(static)
      for (IndexType iq = 0; iq < GetDofHandler().nelements(); ++iq) {
        Transformation(T, iq);
        finiteelement.ReInit(T);

        GlobalToLocal(__U, u, iq);
        GlobalToLocalData(iq, __QN, __QC);

        EQCP->point_cell(GetDofHandler().material(iq));
        integrator.Matrix(*EQCP, __E, finiteelement, __U, __QN, __QC);
        LocalToGlobal(A, __E, iq, d);
      }
      delete EQCP;
      // HANGING NODES
      //    HN->MatrixDiag(u.ncomp(),A);
    }
  }

  ////////////////////////////////////////////////// Integration on the Boundary
  void BoundaryForm(GlobalVector& f,
                    const GlobalVector& u,
                    const ProblemDescriptorInterface& PD,
                    double d) const
  {
    // Do we have a boundary equation?
    if (PD.GetBoundaryEquation() == NULL)
      return;

    LocalParameterData QP;
    GlobalToGlobalData(QP);

    nmatrix<double> T;
    FINITEELEMENT finiteelement;
    INTEGRATOR integrator;
    integrator.BasicInit();
    LocalVector __U, __F;
    LocalData __QN, __QC;

#pragma omp parallel private(T, finiteelement, integrator, __U, __F, __QN, __QC)
    {
      const auto* BEQ = PD.GetBoundaryEquation()->createNew();
      auto COLS = PD.GetBoundaryManager()->GetBoundaryEquationColors();
      BEQ->SetParameterData(QP);

      for (auto col : COLS) {
        const IndexVector& q = *GetDofHandler().ElementOnBoundary(col);
        const IndexVector& l = *GetDofHandler().ElementLocalOnBoundary(col);

#pragma omp for
        for (IndexType i = 0; i < q.size(); i++) {
          IndexType iq = q[i];
          IndexType ile = l[i];

          Transformation(T, iq);
          finiteelement.ReInit(T);

          GlobalToLocal(__U, u, iq);
          GlobalToLocalData(iq, __QN, __QC);

          integrator.BoundaryForm(
            *BEQ, __F, finiteelement, __U, ile, col, __QN, __QC);
          LocalToGlobal(f, __F, iq, d);
        }
      }
      delete BEQ;
    }
  }
  void BoundaryMatrix(MatrixInterface& A,
                      const GlobalVector& u,
                      const ProblemDescriptorInterface& PD,
                      double d) const
  {
    // Do we have a boundary equation?
    if (PD.GetBoundaryEquation() == NULL)
      return;

    LocalParameterData QP;
    GlobalToGlobalData(QP);

    nmatrix<double> T;
    FINITEELEMENT finiteelement;
    INTEGRATOR integrator;
    integrator.BasicInit();
    LocalVector __U;
    LocalData __QN, __QC;
    EntryMatrix __E;

#pragma omp parallel private(T, finiteelement, integrator, __U, __QN, __QC, __E)
    {
      auto* BEQ = PD.GetBoundaryEquation()->createNew();
      const auto COLS = PD.GetBoundaryManager()->GetBoundaryEquationColors();
      BEQ->SetParameterData(QP);

      for (const auto col : COLS) {
        const IndexVector& q = *GetDofHandler().ElementOnBoundary(col);
        const IndexVector& l = *GetDofHandler().ElementLocalOnBoundary(col);
#pragma omp for
        for (IndexType i = 0; i < q.size(); i++) {
          IndexType iq = q[i];
          IndexType ile = l[i];

          Transformation(T, iq);
          finiteelement.ReInit(T);

          GlobalToLocal(__U, u, iq);
          GlobalToLocalData(iq, __QN, __QC);

          integrator.BoundaryMatrix(
            *BEQ, __E, finiteelement, __U, ile, col, __QN, __QC);
          LocalToGlobal(A, __E, iq, d);
        }
      }
      delete BEQ;
    }
  }
  double ComputeBoundaryFunctional(const GlobalVector& u,
                                   const IntSet& Colors,
                                   const BoundaryFunctional& BF) const
  {
    LocalParameterData QP;
    GlobalToGlobalData(QP);

    nmatrix<double> T;

    FINITEELEMENT finiteelement;
    INTEGRATOR integrator;
    integrator.BasicInit();
    LocalVector __U;
    LocalData __QN, __QC;

    double j = 0.;
    for (const auto col : Colors) {
      const IndexVector& q = *GetDofHandler().ElementOnBoundary(col);
      const IndexVector& l = *GetDofHandler().ElementLocalOnBoundary(col);
      for (IndexType i = 0; i < q.size(); i++) {
        IndexType iq = q[i];
        IndexType ile = l[i];

        Transformation(T, iq);
        finiteelement.ReInit(T);

        GlobalToLocal(__U, u, iq);
        GlobalToLocalData(iq, __QN, __QC);

        j += integrator.ComputeBoundaryFunctional(
          BF, finiteelement, ile, col, __U, __QN, __QC);
      }
    }

    return j;
  }

  ////////////////////////////////////////////////// Functionals
  double ComputePointValue(const GlobalVector& u,
                           const Vertex2d& p0,
                           IndexType comp) const
  {
    assert(DIM == 2);
    // very simple version. Only finds nodes
    for (IndexType n = 0; n < GetDofHandler().nnodes(); ++n) {
      double dist = 0;
      for (IndexType d = 0; d < DIM; ++d)
        dist += pow(GetDofHandler().vertex(n)[d] - p0[d], 2.0);
      if (dist < sqrt(1.e-13))
        return u(n, comp);
    }
    std::cerr << "DofHandler::ComputePointValue. Vertex " << p0 << " not found!"
              << std::endl;
    abort();
  }
  double ComputePointValue(const GlobalVector& u,
                           const Vertex3d& p0,
                           IndexType comp) const
  {
    assert(DIM == 2);
    // very simple version. Only finds nodes
    for (IndexType n = 0; n < GetDofHandler().nnodes(); ++n) {
      double dist = 0;
      for (IndexType d = 0; d < DIM; ++d)
        dist += pow(GetDofHandler().vertex(n)[d] - p0[d], 2.0);
      if (dist < sqrt(1.e-13))
        return u(n, comp);
    }
    std::cerr << "DofHandler::ComputePointValue. Vertex " << p0 << " not found!"
              << std::endl;
    abort();
  }
  double ComputePointFunctional(const GlobalVector& u,
                                const PointFunctional& FP) const
  {
    LocalParameterData QP;
    GlobalToGlobalData(QP);
    FP.SetParameterData(QP);

    IndexType dim = GetDofHandler().dimension();
    assert(dim == DIM);
    std::vector<ShortIndexType> comps = FP.GetComps();
    IndexType nn = comps.size();

    std::vector<double> up(nn, 0);

    if (dim == 2) {
      auto v2d = FP.GetPoints2d();
      assert(nn == v2d.size());

      for (IndexType i = 0; i < nn; ++i) {
        up[i] = ComputePointValue(u, v2d[i], comps[i]);
      }
    } else if (dim == 3) {
      auto v3d = FP.GetPoints3d();
      assert(nn == v3d.size());
      for (IndexType i = 0; i < nn; ++i) {
        up[i] = ComputePointValue(u, v3d[i], comps[i]);
      }
    } else {
      std::cout << "wronng dimension: dim = " << dim << std::endl;
      abort();
    }

    return FP.J(up);
  }
  double ComputeDomainFunctional(const GlobalVector& u,
                                 const DomainFunctional& F) const
  {
    LocalParameterData QP;
    GlobalToGlobalData(QP);
    F.SetParameterData(QP);

    nmatrix<double> T;
    FINITEELEMENT finiteelement;
    INTEGRATOR integrator;
    integrator.BasicInit();
    LocalVector __U;
    LocalData __QN, __QC;
    double j = 0.;
    for (IndexType iq = 0; iq < GetDofHandler().nelements(); ++iq) {
      Transformation(T, iq);
      finiteelement.ReInit(T);

      GlobalToLocal(__U, u, iq);
      GlobalToLocalData(iq, __QN, __QC);
      F.point_cell(GetDofHandler().material(iq));
      j +=
        integrator.ComputeDomainFunctional(F, finiteelement, __U, __QN, __QC);
    }
    return j;
  }

  ////////////////////////////////////////////////// Pressure Filter, set
  /// average zero
  void InitFilter(nvector<double>& F) const
  {
    PressureFilter* PF = static_cast<PressureFilter*>(&F);
    (void)PF;
    assert(PF);
    assert(!PF->Active());
  }

  ////////////////////////////////////////////////// Dirichlet Data
  //// NEW Interface
  ////////////////////////////////////////////////// Dirichlet Data
  void StrongDirichletVector(GlobalVector& u,
                             const DirichletData* DD,
                             double d) const
  {
    if (DD == NULL)
      return; // no Dirichlet Data

    for (const auto color : DD->dirichlet_colors()) {
      const IndexVector& bv = *(GetDofHandler().VertexOnBoundary(color));

      LocalParameterData QP;
      GlobalToGlobalData(QP);
      DD->SetParameterData(QP);

      FemData QH;
      nvector<double> ff(u.ncomp(), 0.);

      // threadsafe???
      // kann Dirichlet-Data nicht sowas wie Add Node zeugs haben?
#pragma omp parallel for private(QH) firstprivate(ff)
      for (IndexType i = 0; i < bv.size(); ++i) {
        IndexType index = bv[i];
        QH.clear();
        auto p = GetDataContainer().GetNodeData().begin();
        for (; p != GetDataContainer().GetNodeData().end(); p++) {
          QH[p->first].resize(p->second->ncomp());
          for (IndexType c = 0; c < p->second->ncomp(); c++) {
            QH[p->first][c].m() = p->second->operator()(index, c);
          }
        }
        const Vertex<DIM>& v = GetDofHandler().vertex(index);
        (*DD)(ff, v, color);
        for (auto comp : DD->components_on_color(color))
          u(index, comp) = d * ff[comp];
      }
    }
  }
  void StrongDirichletMatrix(MatrixInterface& A,
                             const ProblemDescriptorInterface& PD) const
  {
    const DirichletData* DD = PD.GetDirichletData();
    if (DD == NULL)
      return; // no Dirichlet Data

    for (auto color : DD->dirichlet_colors()) {
      const IndexVector& bv = *GetDofHandler().VertexOnBoundary(color);
#pragma omp parallel for
      for (IndexType i = 0; i < bv.size(); i++)
        A.dirichlet(bv[i], DD->components_on_color(color));
    }
  }
  void StrongDirichletMatrixOnlyRow(MatrixInterface& A,
                                    const ProblemDescriptorInterface& PD) const
  {
    const DirichletData* DD = PD.GetDirichletData();
    if (DD == NULL)
      return; // no Dirichlet Data

    for (auto color : DD->dirichlet_colors()) {
      const IndexVector& bv = *GetDofHandler().VertexOnBoundary(color);
#pragma omp parallel for
      for (IndexType i = 0; i < bv.size(); i++)
        A.dirichlet_only_row(bv[i], DD->components_on_color(color));
    }
  }
  void StrongDirichletVectorZero(GlobalVector& u,
                                 const ProblemDescriptorInterface& PD) const
  {
    const DirichletData* DD = PD.GetDirichletData();
    if (DD == NULL)
      return; // no Dirichlet Data

    for (auto color : DD->dirichlet_colors()) {
      const IndexVector& bv = *GetDofHandler().VertexOnBoundary(color);

#pragma omp parallel for
      for (IndexType i = 0; i < bv.size(); i++) {
        IndexType index = bv[i];
        for (auto comp : DD->components_on_color(color))
          u(index, comp) = 0.;
      }
    }
  }

  ////////////////////////////////////////////////// old interface
  ////////////////////////////////////////////////// Dirichlet Data
  void StrongDirichletVector(GlobalVector& u,
                             const DirichletData& BF,
                             IndexType col,
                             const std::vector<IndexType>& comp,
                             double d) const
  {
    abort();
  }
  void StrongDirichletMatrix(MatrixInterface& A,
                             IndexType col,
                             const std::vector<IndexType>& comp) const
  {
    abort();
  }
  void StrongDirichletMatrixOnlyRow(MatrixInterface& A,
                                    IndexType col,
                                    const std::vector<IndexType>& comp) const
  {
    abort();
  }

  ////////////////////////////////////////////////// Errors
  void ComputeError(const GlobalVector& u,
                    LocalVector& err,
                    const ExactSolution* ES) const
  {
    IndexType ncomp = u.ncomp();
    err.ncomp() = ncomp;
    err.reservesize(3);
    err = 0.;

    GlobalVector lerr(ncomp, 3);
    lerr.zero();

    LocalParameterData QP;
    GlobalToGlobalData(QP);
    ES->SetParameterData(QP);

    nmatrix<double> T;
    LocalVector U;
    LocalData QN, QC;
    FINITEELEMENT finiteelement;
    INTEGRATOR integrator;
    integrator.BasicInit();
    for (IndexType iq = 0; iq < GetDofHandler().nelements(); iq++) {
      Transformation(T, iq);
      finiteelement.ReInit(T);
      GlobalToLocal(U, u, iq);
      GlobalToLocalData(iq, QN, QC);
      integrator.ErrorsByExactSolution(lerr, finiteelement, *ES, U, QN, QC);

      for (IndexType c = 0; c < ncomp; c++) {
        err(0, c) += lerr(0, c);
        err(1, c) += lerr(1, c);
        err(2, c) = std::max(err(2, c), lerr(2, c));
      }
    }
    for (IndexType c = 0; c < ncomp; c++) {
      err(0, c) = sqrt(err(0, c));
      err(1, c) = sqrt(err(1, c));
    }
  }
};

// #include "finiteelement.h"
// #include "elementintegrator.h"
// #include "integrationformula.h"
// #include "transformation2d.h"
// #include "transformation3d.h"
// #include "baseq22d.h"
// #include "baseq23d.h"

// namespace Gascoigne
// {

// QP Finite Elements
#define CGBaseQ12d CGBase<2, 2>
#define CGBaseQ22d CGBase<2, 3>
#define CGBaseQ42d CGBase<2, 5>

#define CGTransformationQ12d Transformation2d<CGBaseQ12d>
#define CGTransformationQ22d Transformation2d<CGBaseQ22d>
#define CGTransformationQ42d Transformation2d<CGBaseQ42d>

#define LagrangeDiscQ12d                                                       \
  LagrangeDisc<2,                                                              \
               2,                                                              \
               FiniteElement<2, 1, CGTransformationQ12d, CGBaseQ12d>,          \
               ElementIntegratorQ12d>
#define LagrangeDiscQ22d                                                       \
  LagrangeDisc<2,                                                              \
               3,                                                              \
               FiniteElement<2, 1, CGTransformationQ22d, CGBaseQ22d>,          \
               ElementIntegratorQ22d>
#define LagrangeDiscQ42d                                                       \
  LagrangeDisc<2,                                                              \
               5,                                                              \
               FiniteElement<2, 1, CGTransformationQ42d, CGBaseQ42d>,          \
               ElementIntegratorQ42d>

#define CGBaseQ13d CGBase<3, 2>
#define CGBaseQ23d CGBase<3, 3>
#define CGBaseQ43d CGBase<3, 5>

#define CGTransformationQ13d Transformation3d<CGBaseQ13d>
#define CGTransformationQ23d Transformation3d<CGBaseQ23d>
#define CGTransformationQ43d Transformation3d<CGBaseQ43d>

#define LagrangeDiscQ13d                                                       \
  LagrangeDisc<3,                                                              \
               2,                                                              \
               FiniteElement<3, 2, CGTransformationQ13d, CGBaseQ13d>,          \
               ElementIntegratorQ13d>
#define LagrangeDiscQ23d                                                       \
  LagrangeDisc<3,                                                              \
               3,                                                              \
               FiniteElement<3, 2, CGTransformationQ23d, CGBaseQ23d>,          \
               ElementIntegratorQ23d>
#define LagrangeDiscQ43d                                                       \
  LagrangeDisc<3,                                                              \
               5,                                                              \
               FiniteElement<3, 2, CGTransformationQ43d, CGBaseQ43d>,          \
               ElementIntegratorQ43d>

} // namespace Gascoigne

/*----------------------------   lagrangedisc.h ---------------------------*/
/* end of #ifndef __lagrangedisc_H */
#endif
/*----------------------------   lagrangedisc.h ---------------------------*/
