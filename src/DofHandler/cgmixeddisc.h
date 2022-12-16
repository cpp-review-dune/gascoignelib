/*----------------------------   cgdisc.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __cgmixeddisc_H
#define __cgmixeddisc_H
/*----------------------------   cgdisc.h     ---------------------------*/

/**
 *
 * Copyright (C) 2018 by the Gascoigne 3D authors
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
//#include "omp.h"

#include "../Common/stopwatch.h"
#include "../Discretization/Q1/baseq12d.h"
#include "../Discretization/Q1/baseq13d.h"
#include "../Discretization/Q1/finiteelement.h"
#include "../Discretization/Q1/gascoignevisualization.h"
#include "../Discretization/Q1/hnstructureq13d.h"
#include "../Discretization/Q1/mginterpolatornested.h"
#include "../Discretization/Q1/pressurefilter.h"
#include "../Discretization/Q2/baseq13dpatch.h"
#include "../Discretization/Q2/baseq1patch.h"
#include "../Discretization/Q2/baseq22d.h"
#include "../Discretization/Q2/baseq23d.h"
#include "../Discretization/Q2/hnstructureq22d.h"
#include "../Discretization/Q2/hnstructureq23d.h"
#include "../Interface/discretizationinterface.h"
#include "../LinAlg/sparsestructure.h"
#include "../Problems/problemdescriptorbase.h"

#include "dofhandler.h"
#include "elementlpsintegrator.h"
#include "mixedelementintegrator.h"

namespace Gascoigne {
namespace atom_ops {
inline void
add_node_mixed(double s,
               IndexType i_f,
               GlobalVector& __restrict__ f,
               IndexType i_F,
               const LocalVector& __restrict__ F)
{
  const IndexType iif = i_f * f.ncomp();
  const IndexType iiF = i_F * F.ncomp();
  for (IndexType c = 0; c < f.ncomp(); c++) {
#pragma omp atomic update
    f[iif + c] += s * F[iiF + c];
  }
}
} // namespace atom_ops

/////////////////////////////////////////////
///
///@brief
///  ... comments DiscretizationInterface

///
///
/////////////////////////////////////////////

// DIM=2,3
// DEGREE = 1 (Q1) 2 (Q2)
template<int DIM,
         int DEGREE,
         class FINITEELEMENTTRIAL,
         class FINITEELEMENTTEST,
         class INTEGRATOR>
class CGMixedDisc : public DiscretizationInterface
{
private:
  const DofHandler<DIM>* dofhandler;
  mutable DataContainer datacontainer;

protected:
  // Hanging nodes
  HNStructureInterface* HN;

public:
  CGMixedDisc()
    : HN(NULL)
  {}
  ~CGMixedDisc() {}

  //    HNStructureInterface* NewHNStructure() {abort();}
  HNStructureInterface* NewHNStructure()
  {
    // This will not work at the moment.
    // We need HNStructureQ2 and HNStructureQ1Patch
    if (DEGREE == 1) {
      if (DIM == 3)
        return new HNStructureQ13d;
      else if (DIM == 2)
        return new HNStructureQ12d;
      else
        assert(0);
    } else if (DEGREE == 2) {
      if (DIM == 3)
        return new HNStructureQ23d;
      else if (DIM == 2)
        return new HNStructureQ22d;
      else
        assert(0);
    } else
      abort();
  }

  const DofHandler<DIM>* GetDofHandler() const { return dofhandler; }

  const DataContainer& GetDataContainer() const { return datacontainer; }
  void SetDataContainer(const DataContainer& q) { datacontainer = q; }

  //
  //// Functions called from the Solver
  //
  std::string GetName() const { return "CG Mixed Discretization"; }

  // Q12d Q13d Q22d Q32d

  // Visualization
  void VisuVtk(const ComponentInformation* CI,
               const ParamFile& pf,
               const std::string& name,
               const GlobalVector& u,
               IndexType i) const
  {
    HNAverage(const_cast<GlobalVector&>(u));

    GascoigneVisualization Visu;
    Visu.SetMesh(GetDofHandler());
    if (CI) {
      Visu.AddPointVector(CI, &u);
    } else {
      Visu.AddPointVector(&u);
    }

    Visu.read_parameters(pf);
    Visu.set_name(name);
    Visu.step(i);
    Visu.write();

    HNZero(const_cast<GlobalVector&>(u));
  }

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
    assert(HN == NULL);
    HN = NewHNStructure();
    assert(HN);
  }
  void ReInit(const GascoigneMesh* M)
  {
    dofhandler = dynamic_cast<const DofHandler<DIM>*>(M);
    assert(dofhandler);

    assert(HN);
    HN->ReInit(M);
  }

  Vertex2d vertex2d(IndexType i) const
  {
    assert(i < GetDofHandler()->nnodes());
    return GetDofHandler()->vertex2d(i);
  }
  Vertex3d vertex3d(IndexType i) const
  {
    assert(i < GetDofHandler()->nnodes());
    return GetDofHandler()->vertex3d(i);
  }

  IndexType ndofs() const { return GetDofHandler()->nnodes(); }
  IndexType nelements() const { return GetDofHandler()->nelements(DEGREE); }
  IndexType ndegree() const { return DEGREE; }
  IndexType nhanging() const { return GetDofHandler()->nhanging(); }

  IndexType ndofs_withouthanging() const { return ndofs() - nhanging(); }

  // Hanging nodes
  void HNAverage(GlobalVector& x) const
  {
    assert(HN);
    HN->Average(x);
  }
  void HNDistribute(GlobalVector& x) const
  {
    assert(HN);
    HN->Distribute(x);
  }
  void HNZero(GlobalVector& x) const
  {
    assert(HN);
    HN->Zero(x);
  }
  bool HNZeroCheck(const GlobalVector& x) const
  {
    assert(HN);
    return HN->ZeroCheck(x);
  }
  void HNAverageData() const
  {
    const GlobalData& gd = GetDataContainer().GetNodeData();
    for (const auto& it : gd)
      HNAverage(*const_cast<GlobalVector*>(it.second));
  }
  void HNZeroData() const
  {
    const GlobalData& gd = GetDataContainer().GetNodeData();
    for (const auto& it : gd)
      HNZero(*const_cast<GlobalVector*>(it.second));
  }

  //////////////////////////////////////////////////
  virtual void Transformation(nmatrix<double>& T, IndexType iq) const
  {
    assert(GetDofHandler()->dimension() == DIM);
    IndexType ne = GetDofHandler()->nodes_per_element(DEGREE);

    IndexVector indices = GetDofHandler()->GetElement(DEGREE, iq);
    assert(ne == indices.size());

    T.memory(DIM, ne);
    if (DIM == 2) {
      for (IndexType ii = 0; ii < ne; ii++) {
        Vertex2d v = GetDofHandler()->vertex2d(indices[ii]);
        T(0, ii) = v.x();
        T(1, ii) = v.y();
      }
    } else if (DIM == 3) {
      for (IndexType ii = 0; ii < ne; ii++) {
        Vertex3d v = GetDofHandler()->vertex3d(indices[ii]);
        T(0, ii) = v.x();
        T(1, ii) = v.y();
        T(2, ii) = v.z();
      }
    }
  }
  void ConstructInterpolator(MgInterpolatorInterface* I,
                             const MeshTransferInterface* MT)
  {
    MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
    assert(IP);
    IP->BasicInit(MT);
  }

  void Structure(SparseStructureInterface* SI) const
  {
    SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
    assert(S);

    S->build_begin(ndofs());
    for (IndexType iq = 0; iq < GetDofHandler()->nelements(DEGREE); iq++) {
      IndexVector indices = GetDofHandler()->GetElement(DEGREE, iq);
      // HANGING NODES
      HN->CondenseHanging(indices);
      S->build_add(indices.begin(), indices.end());
    }
    S->build_end();
    // HANGING NODES
    HN->SparseStructureDiag(S);
  }

  ////////////////////////////////////////////////// handling local / global
  virtual void GlobalToGlobalData(LocalParameterData& QP) const
  {
    const GlobalParameterData& gpd = GetDataContainer().GetParameterData();
    QP.clear();

    for (auto p : gpd)
      QP.insert(make_pair(p.first, *p.second));
  }
  // virtual void LocalToGlobal_ohnecritic(MatrixInterface& A,
  //                                       EntryMatrix& E,
  //                                       IndexType iq,
  //                                       double s) const
  // {
  //   IntVector indices = GetDofHandler()->GetElement(DEGREE, iq);

  //   // HANGING NODES
  //   HN->CondenseHanging(E, indices);
  //   IntVector::const_iterator start = indices.begin();
  //   IntVector::const_iterator stop = indices.end();
  //   //#pragma omp critical
  //   A.entry(start, stop, E, s);
  // }
  // virtual void LocalToGlobal_ohnecritic(GlobalVector& f,
  //                                       const LocalVector& F,
  //                                       IndexType iq,
  //                                       double s) const
  // {
  //   IntVector indices = GetDofHandler()->GetElement(DEGREE, iq);
  //   for (IndexType ii = 0; ii < indices.size(); ii++) {
  //     IndexType i = indices[ii];
  //     //#pragma omp critical
  //     f.add_node_mixed(i, s, ii, F);
  //   }
  // }

  virtual void LocalToGlobal(MatrixInterface& A,
                             EntryMatrix& E,
                             IndexType iq,
                             double s) const
  {
    IndexVector indices = GetDofHandler()->GetElement(DEGREE, iq);

    // HANGING NODES
    HN->CondenseHanging(E, indices);
    IndexVector::const_iterator start = indices.begin();
    IndexVector::const_iterator stop = indices.end();

    A.entry(start, stop, E, s);
  }
  virtual void LocalToGlobal(GlobalVector& f,
                             const LocalVector& F,
                             IndexType iq,
                             double s) const
  {
    IndexVector indices = GetDofHandler()->GetElement(DEGREE, iq);
    for (IndexType ii = 0; ii < indices.size(); ii++) {
      atom_ops::add_node_mixed(s, indices[ii], f, ii, F);
    }
  }
  virtual void GlobalToLocal(LocalVector& U,
                             const GlobalVector& u,
                             IndexType iq) const
  {
    IndexVector indices = GetDofHandler()->GetElement(DEGREE, iq);
    U.ReInit(u.ncomp(), indices.size());
    for (IndexType ii = 0; ii < indices.size(); ii++) {
      IndexType i = indices[ii];
      U.equ_node(ii, i, u);
    }
  }
  virtual void GlobalToLocalCell(LocalVector& U,
                                 const GlobalVector& u,
                                 IndexType iq) const
  {
    U.ReInit(u.ncomp(), 1);
    for (IndexType c = 0; c < u.ncomp(); ++c)
      U(0, c) = u(iq, c);
  }
  virtual void GlobalToLocalData(IndexType iq,
                                 LocalData& QN,
                                 LocalData& QC) const
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
    //!!!
    u.zero();
  }

  // assemble of the weak formulation for all test functions
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
      FINITEELEMENTTRIAL finiteelementtrial;
      FINITEELEMENTTEST finiteelementtest;

      INTEGRATOR integrator;
      integrator.BasicInit();

      LocalVector __U, __F;
      LocalData __QN, __QC;
      const auto EQCP = EQ.createNew();

      EQCP->SetParameterData(QP);
      // EQCP->SetTime(PD.time(),PD.dt());

#pragma omp for schedule(static)
      for (IndexType iq = 0; iq < GetDofHandler()->nelements(DEGREE); ++iq) {
        Transformation(T, iq);
        finiteelementtrial.ReInit(T);
        finiteelementtest.ReInit(T);

        GlobalToLocal(__U, u, iq);
        GlobalToLocalData(iq, __QN, __QC);

        EQCP->point_cell(GetDofHandler()->material(DEGREE, iq));

        integrator.Form(
          *EQCP, __F, finiteelementtrial, finiteelementtest, __U, __QN, __QC);

        LocalToGlobal(f, __F, iq, d);
      }
      delete EQCP;
    }
  }

  void DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const
  {
    abort();
  }

  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
  {
    LocalParameterData QP;
    GlobalToGlobalData(QP);
#pragma omp parallel
    {
      nmatrix<double> T;
      FINITEELEMENTTRIAL finiteelementtrial;
      FINITEELEMENTTEST finiteelementtest;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __F;
      LocalData __QN, __QC;
      const auto RHSCP = RHS.createNew();
      RHSCP->SetParameterData(QP);
      // RHSCP->SetTime(PD.time(),PD.dt());

#pragma omp for schedule(static)
      for (IndexType iq = 0; iq < GetDofHandler()->nelements(DEGREE); ++iq) {
        Transformation(T, iq);

        finiteelementtrial.ReInit(T);
        finiteelementtest.ReInit(T);

        GlobalToLocalData(iq, __QN, __QC);
        RHSCP->point_cell(GetDofHandler()->material(DEGREE, iq));
        integrator.Rhs(
          *RHSCP, __F, finiteelementtrial, finiteelementtest, __QN, __QC);
        LocalToGlobal(f, __F, iq, s);
      }
      delete RHSCP;
    }
  }

  void BoundaryRhs(GlobalVector& f,
                   const IntSet& Colors,
                   const BoundaryRightHandSide& BRHS,
                   double s) const
  {
    abort();
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
      FINITEELEMENTTRIAL finiteelementtrial;
      FINITEELEMENTTEST finiteelementtest;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __U;
      LocalData __QN, __QC;
      EntryMatrix __E;

      const auto EQCP = EQ.createNew();
      EQCP->SetParameterData(QP);
      // EQ->SetTime(PD.time(),PD.dt());
#pragma omp for schedule(static)
      for (IndexType iq = 0; iq < GetDofHandler()->nelements(DEGREE); ++iq) {
        Transformation(T, iq);
        finiteelementtrial.ReInit(T);
        finiteelementtest.ReInit(T);

        GlobalToLocal(__U, u, iq);
        GlobalToLocalData(iq, __QN, __QC);

        EQCP->point_cell(GetDofHandler()->material(DEGREE, iq));
        // EQ.cell(GetDofHandler(),iq,__U,__QN);
        integrator.Matrix(
          *EQCP, __E, finiteelementtrial, finiteelementtest, __U, __QN, __QC);
        LocalToGlobal(A, __E, iq, d);
      }
      delete EQCP;
    }
    HN->MatrixDiag(u.ncomp(), A);
  }

  ////////////////////////////////////////////////// Integration on the Boundary
  void BoundaryForm(GlobalVector& f,
                    const GlobalVector& u,
                    const ProblemDescriptorInterface& PD,
                    double d) const
  {}
  void BoundaryMatrix(MatrixInterface& A,
                      const GlobalVector& u,
                      const ProblemDescriptorInterface& PD,
                      double d) const
  {}

  ////////////////////////////////////////////////// Functionals

  // Computes the divergence (squared) on one element given by Vector M
  double LocalDiv(const LocalVector& U, const LocalVector& M) const { abort(); }
  double ComputeBoundaryFunctional(const GlobalVector& u,
                                   const IntSet& Colors,
                                   const BoundaryFunctional& BF) const
  {
    abort();
  }

  void VertexTransformation(const Vertex<DIM>& p0,
                            Vertex<DIM>& p,
                            IndexType iq) const
  {
    abort();
  }
  IndexType GetElementNumber(const Vertex<DIM>& p0, Vertex<DIM>& p) const
  {
    abort();
  }

  ////////////////////////////////////////////////// Functionals
  double ComputePointValue(const GlobalVector& u,
                           const Vertex2d& p0,
                           IndexType comp) const
  {
    abort();
  }
  double ComputePointValue(const GlobalVector& u,
                           const Vertex3d& p0,
                           IndexType comp) const
  {
    abort();
  }
  double ComputePointFunctional(const GlobalVector& u,
                                const PointFunctional& FP) const
  {
    abort();
  }

  double ComputeDomainFunctional(const GlobalVector& u,
                                 const DomainFunctional& F) const
  {
    abort();
  }
  void MassMatrix(MatrixInterface& A) const { abort(); }

  void InitFilter(nvector<double>& F) const {}
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
      const IndexVector& bv = *GetDofHandler()->VertexOnBoundary(color);

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
        const Vertex<DIM>& v = GetDofHandler()->vertex(index);
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
      const IndexVector& bv = *GetDofHandler()->VertexOnBoundary(color);
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
      const IndexVector& bv = *GetDofHandler()->VertexOnBoundary(color);
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
      const IndexVector& bv = *GetDofHandler()->VertexOnBoundary(color);

#pragma omp parallel for
      for (IndexType i = 0; i < bv.size(); i++) {
        IndexType index = bv[i];
        for (auto comp : DD->components_on_color(color))
          u(index, comp) = 0.;
      }
    }
  }

  void StrongPeriodicVector(GlobalVector& u,
                            const PeriodicData& BF,
                            IndexType col,
                            const std::vector<IndexType>& comp,
                            double d) const
  {
    const GascoigneMesh* GMP =
      dynamic_cast<const GascoigneMesh*>(GetDofHandler());
    assert(GMP);
    DoubleVector ff(u.ncomp(), 0.);
    const IndexVector& bv = *GMP->VertexOnBoundary(col);

    FemData QH;

    LocalParameterData QP;
    GlobalToGlobalData(QP);
    BF.SetParameterData(QP);

    // for(IndexType ii=0;ii<comp.size();ii++)
    //{
    //  IndexType c = comp[ii];
    //  if(c<0) {
    //    cerr << "negative component: " << c << endl;
    //    abort();
    //  } else if(c>=u.ncomp()){
    //    cerr << "unknown component: " << c << endl;
    //    abort();
    //  }
    //}

    for (IndexType i = 0; i < bv.size(); i++) {
      IndexType index = bv[i];

      QH.clear();
      GlobalData::const_iterator p = GetDataContainer().GetNodeData().begin();
      for (; p != GetDataContainer().GetNodeData().end(); p++) {
        QH[p->first].resize(p->second->ncomp());
        for (IndexType c = 0; c < p->second->ncomp(); c++) {
          QH[p->first][c].m() = p->second->operator()(index, c);
        }
      }

      BF.SetFemData(QH);

      const Vertex2d& v = GMP->vertex2d(index);

      BF(ff, v, col);
      for (IndexType iii = 0; iii < comp.size(); iii++) {
        IndexType c = comp[iii];
        u(index, c) = d * ff[c];
      }
    }
  }

  ////////////////////////////////////////////////// Errors
  void ComputeError(const GlobalVector& u,
                    LocalVector& err,
                    const ExactSolution* ES) const
  {
    //    abort();
  }
};

////// LPS
typedef CGMixedDisc<
  2,
  2,
  FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>, BaseQ12dPatch>,
  FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>, BaseQ22d>,
  MixedElementIntegratorQ12dPatch>
  CGMixedDiscQ12dPatch;

} // namespace Gascoigne

/*----------------------------   cgdisc.h     ---------------------------*/
/* end of #ifndef __cgdisc_H */
#endif
/*----------------------------   cgdisc.h     ---------------------------*/
