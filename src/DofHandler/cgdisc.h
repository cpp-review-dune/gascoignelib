/*----------------------------   cgdisc.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __cgdisc_H
#define __cgdisc_H
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

#include "discretizationinterface.h"
#include "mginterpolatornested.h"
//#include "omp.h"
#include "pressurefilter.h"
#include "problemdescriptorbase.h"
#include "sparsestructure.h"
#include "stopwatch.h"
#include "gascoignevisualization.h"

#include "hnstructureq22d.h"
#include "hnstructureq23d.h"
#include "hnstructureq22d.h"


namespace Gascoigne
{

  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments DiscretizationInterface

  ///
  ///
  /////////////////////////////////////////////

  // DIM=2,3
  // DEGREE = 1 (Q1) 2 (Q2)
  template <int DIM, int DEGREE, class FINITEELEMENT, class INTEGRATOR>
  class CGDisc : public DiscretizationInterface
  {
  private:
    const DofHandler<DIM> *dofhandler;
    mutable DataContainer datacontainer;

  protected:


    // Hanging nodes
    HNStructureInterface*    HN;
    std::vector<std::vector<int> >  Coloring;	
    
    
  public:
    CGDisc():HN(NULL) {}
    ~CGDisc() {}

    //    HNStructureInterface* NewHNStructure() {abort();}
    HNStructureInterface* NewHNStructure()
    {
      if (DIM==3)
	return new HNStructureQ23d;
      else if (DIM==2)
	return new HNStructureQ22d;
      else assert(0);
    }

    const DofHandler<DIM> *GetDofHandler() const
    {
      return dofhandler;
    }


    const DataContainer &GetDataContainer() const
    {
      return datacontainer;
    }
    void SetDataContainer(const DataContainer &q)
    {
      datacontainer = q;
    }


    //
    //// Functions called from the Solver
    //
    std::string GetName() const
    {
      return "CG Discretization";
    }
    // Visualization
    void VisuVtk(const ComponentInformation* CI, const ParamFile& pf,
		 const std::string &name, const GlobalVector& u, int i) const
    {
      HNAverage(const_cast<GlobalVector &>(u));
      
      GascoigneVisualization Visu;
      
      Visu.SetMesh(GetDofHandler());
      
      if (CI)
	{
	  Visu.AddPointVector(CI, &u);
	}
      else
      {
	Visu.AddPointVector(&u);
      }
      
      Visu.read_parameters(&pf);
      Visu.set_name(name);
      Visu.step(i);
      Visu.write();
      
      HNZero(const_cast<GlobalVector &>(u));
    }

    
    void AddNodeVector(const std::string &name, const GlobalVector *q) const
    {
      datacontainer.AddNodeVector(name, q);
    }
    void DeleteNodeVector(const std::string &name) const
    {
      datacontainer.DeleteNodeVector(name);
    }

    void AddCellVector(const std::string &name, const GlobalVector *q) const
    {
      datacontainer.AddCellVector(name, q);
    }
    void DeleteCellVector(const std::string &name) const
    {
      datacontainer.DeleteCellVector(name);
    }

    void AddParameterVector(const std::string &name,
                            const GlobalParameterVector *q) const
    {
      datacontainer.AddParameterVector(name, q);
    }
    void DeleteParameterVector(const std::string &name) const
    {
      datacontainer.DeleteParameterVector(name);
    }

    void BasicInit(const ParamFile *pf)
    {
      assert(HN==NULL);
      HN = NewHNStructure();
      assert(HN);
    }
    void ReInit(const GascoigneMesh *M)
    {
      dofhandler = dynamic_cast<const DofHandler<DIM> *>(M);
      assert(dofhandler);
		
      assert(HN);
      HN->ReInit(M);
      
	  std::cout<<"huhhhhhuuuuu coloring for dofhandler"<<std::endl;	    
      ElementColoring(DEGREE);
      std::cout<<"GetDofHandler()->NumberofColors(): "<<NumberofColors()<<std::endl;
    }

    
    Vertex2d vertex2d(int i) const{ assert(i<GetDofHandler()->nnodes()); return GetDofHandler()->vertex2d(i); }
    Vertex3d vertex3d(int i) const{ assert(i<GetDofHandler()->nnodes()); return GetDofHandler()->vertex3d(i); }
    

    

    int ndofs() const
    {
      return GetDofHandler()->nnodes();
    }
    int nelements() const
    {
      return GetDofHandler()->nelements(DEGREE);
    }
    int ndofs_withouthanging() const
    {
      return ndofs();
      // HANGING NODES
    }

    // Hanging nodes
    void HNAverage   (GlobalVector& x) const { assert(HN); HN->Average(x); }
    void HNDistribute(GlobalVector& x) const { assert(HN); HN->Distribute(x); }
    void HNZero      (GlobalVector& x) const { assert(HN); HN->Zero(x); }
    bool HNZeroCheck (const GlobalVector& x) const { assert(HN); return HN->ZeroCheck(x); }
    void HNAverageData() const
    {
      const GlobalData& gd = GetDataContainer().GetNodeData();
      for (const auto &it : gd)
	HNAverage(*const_cast<GlobalVector*>(it.second));
    }
    void HNZeroData   () const
    {
      const GlobalData& gd = GetDataContainer().GetNodeData();
      for (const auto &it : gd)
	HNZero(*const_cast<GlobalVector*>(it.second));
    }


    //////////////////////////////////////////////////
    void Transformation(nmatrix<double> &T, int iq) const
    {
      assert(GetDofHandler()->dimension() == DIM);
      int ne = GetDofHandler()->nodes_per_element(DEGREE);

      IntVector indices = GetDofHandler()->GetElement(DEGREE, iq);
      assert(ne == indices.size());

      T.memory(DIM, ne);
      if (DIM == 2)
      {
        for (int ii = 0; ii < ne; ii++)
        {
          Vertex2d v = GetDofHandler()->vertex2d(indices[ii]);
          T(0, ii) = v.x();
          T(1, ii) = v.y();
        }
      }
      else if (DIM == 3)
      {
        for (int ii = 0; ii < ne; ii++)
        {
          Vertex3d v = GetDofHandler()->vertex3d(indices[ii]);
          T(0, ii) = v.x();
          T(1, ii) = v.y();
          T(2, ii) = v.z();
        }
      }
    }
    void ConstructInterpolator(MgInterpolatorInterface *I,
                               const MeshTransferInterface *MT)
    {
      MgInterpolatorNested *IP = dynamic_cast<MgInterpolatorNested *>(I);
      assert(IP);
      IP->BasicInit(MT);
    }


    void Structure(SparseStructureInterface *SI) const
    {
      SparseStructure *S = dynamic_cast<SparseStructure *>(SI);
      assert(S);

      S->build_begin(ndofs());
      for (int iq = 0; iq < GetDofHandler()->nelements(DEGREE); iq++)
      {
        IntVector indices = GetDofHandler()->GetElement(DEGREE, iq);
        // HANGING NODES
	HN->CondenseHanging(indices);
        S->build_add(indices.begin(), indices.end());
      }
      S->build_end();
      // HANGING NODES
      HN->SparseStructureDiag(S);
    }


    ////////////////////////////////////////////////// handling local / global
    void GlobalToGlobalData(LocalParameterData &QP) const
    {
      const GlobalParameterData &gpd = GetDataContainer().GetParameterData();
      QP.clear();

      for (auto p : gpd)
        QP.insert(make_pair(p.first, *p.second));
    }
    void
    LocalToGlobal(MatrixInterface &A, EntryMatrix &E, int iq, double s) const
    {
      IntVector indices = GetDofHandler()->GetElement(DEGREE, iq);

      // HANGING NODES
       HN->CondenseHanging(E,indices);
      IntVector::const_iterator start = indices.begin();
      IntVector::const_iterator stop = indices.end();
#pragma omp critical
      A.entry(start, stop, E, s);
    }
    void
    LocalToGlobal(GlobalVector &f, const LocalVector &F, int iq, double s) const
    {
      IntVector indices = GetDofHandler()->GetElement(DEGREE, iq);
      for (int ii = 0; ii < indices.size(); ii++)
      {
        int i = indices[ii];
#pragma omp critical
        f.add_node(i, s, ii, F);
      }
    }
    void GlobalToLocal(LocalVector &U, const GlobalVector &u, int iq) const
    {
      IntVector indices = GetDofHandler()->GetElement(DEGREE, iq);
      U.ReInit(u.ncomp(), indices.size());
      for (int ii = 0; ii < indices.size(); ii++)
      {
        int i = indices[ii];
        U.equ_node(ii, i, u);
      }
    }
    void GlobalToLocalCell(LocalVector &U, const GlobalVector &u, int iq) const
    {
      U.ReInit(u.ncomp(), 1);
      for (int c = 0; c < u.ncomp(); ++c)
        U(0, c) = u(iq, c);
    }
    void GlobalToLocalData(int iq, LocalData &QN, LocalData &QC) const
    {
      const GlobalData &gnd = GetDataContainer().GetNodeData();
      QN.clear();
      GlobalData::const_iterator p = gnd.begin();
      for (; p != gnd.end(); p++)
      {
        GlobalToLocal(QN[p->first], *p->second, iq);
      }

      const GlobalData &gcd = GetDataContainer().GetCellData();
      QC.clear();
      GlobalData::const_iterator q = gcd.begin();
      for (; q != gcd.end(); q++)
      {
        GlobalToLocalCell(QC[q->first], *q->second, iq);
      }
    }
    void InterpolateSolution(GlobalVector &u, const GlobalVector &uold) const;

    ////////////////////////////////////////////////// integration ueber Zellen   
    void Form(GlobalVector &f,
              const GlobalVector &u,
              const ProblemDescriptorInterface &PD,
              double d) const
    {
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
        auto *EQ = PD.NewEquation();
        EQ->SetParameterData(QP);
        //for(int col=0;col<NumberofColors();col++)
        {
        //const std::vector<int>& ewcol= elementswithcolor(col);
#pragma omp for
		    for (int iq = 0; iq < GetDofHandler()->nelements(DEGREE); ++iq)
		    //for (int iii = 0; iii < ewcol.size(); ++iii)
		    {
		    //  int iq=ewcol[iii];
		      Transformation(T, iq);
		      finiteelement.ReInit(T);

		      GlobalToLocal(__U, u, iq);
		      GlobalToLocalData(iq, __QN, __QC);
		      EQ->point_cell(GetDofHandler()->material(DEGREE, iq));

		      integrator.Form(*EQ, __F, finiteelement, __U, __QN, __QC);
		      LocalToGlobal(f, __F, iq, d);
		    }
        }
        delete EQ;
      }
    }
    void Rhs(GlobalVector &f, const DomainRightHandSide &RHS, double s) const
    {
      LocalParameterData QP;
      GlobalToGlobalData(QP);
      RHS.SetParameterData(QP);


      nmatrix<double> T;
      FINITEELEMENT finiteelement;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __F;
      LocalData __QN, __QC;

#pragma omp parallel for private(T, finiteelement, integrator, __F, __QN, __QC)
      for (int iq = 0; iq < GetDofHandler()->nelements(DEGREE); ++iq)
      {
        Transformation(T, iq);
        finiteelement.ReInit(T);

        GlobalToLocalData(iq, __QN, __QC);
        RHS.point_cell(GetDofHandler()->material(DEGREE, iq));
        integrator.Rhs(RHS, __F, finiteelement, __QN, __QC);
        LocalToGlobal(f, __F, iq, s);
      }
    }
    void Matrix(MatrixInterface &A,
                const GlobalVector &u,
                const ProblemDescriptorInterface &PD,
                double d) const
    {
      LocalParameterData QP;
      GlobalToGlobalData(QP);

      nmatrix<double> T;
      FINITEELEMENT finiteelement;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __U;
      LocalData __QN, __QC;
      EntryMatrix __E;

#pragma omp parallel private(                                                  \
    T, finiteelement, integrator, __U,  __QN, __QC, __E)
      {
        Equation *EQ = PD.NewEquation();
        EQ->SetParameterData(QP);
        //for(int col=0;col<NumberofColors();col++)
        {
        //const std::vector<int>& ewcol= elementswithcolor(col);
#pragma omp for
		    for (int iq = 0; iq < GetDofHandler()->nelements(DEGREE); ++iq)
		    //for (int iii = 0; iii < ewcol.size(); ++iii)
		    {
		      //int iq=ewcol[iii];
        //{
//#pragma omp for
//          for (int iq = 0; iq < GetDofHandler()->nelements(DEGREE); ++iq)
//          {
            Transformation(T, iq);
            finiteelement.ReInit(T);

            GlobalToLocal(__U, u, iq);
            GlobalToLocalData(iq, __QN, __QC);

            EQ->point_cell(GetDofHandler()->material(DEGREE, iq));
            // EQ.cell(GetDofHandler(),iq,__U,__QN);
            integrator.Matrix(*EQ, __E, finiteelement, __U, __QN, __QC);
            LocalToGlobal(A, __E, iq, d);
          }
         }
          // HANGING NODES
	  HN->MatrixDiag(u.ncomp(),A);
          delete EQ;
        //}
      }
    }


    ////////////////////////////////////////////////// Integration on the Boundary
    void BoundaryForm(GlobalVector& f, const GlobalVector& u, const ProblemDescriptorInterface& PD, double d) const 
    {
      // Do we have a boundary equation?
      if (PD.NewBoundaryEquation() == NULL)
	return;
      
      LocalParameterData QP;
      GlobalToGlobalData(QP);
      
      nmatrix<double> T;
      FINITEELEMENT finiteelement;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __U, __F;
      LocalData __QN, __QC;
      
#pragma omp parallel private(T,finiteelement,integrator, __U,__F,__QN,__QC)
      {
	const auto *BEQ = PD.NewBoundaryEquation();
	auto COLS = PD.GetBoundaryManager()->GetBoundaryEquationColors();
	BEQ->SetParameterData(QP);
	
	for(auto col : COLS)
	  {
	    const IntVector& q = *GetDofHandler()->ElementOnBoundary(DEGREE,col);
	    const IntVector& l = *GetDofHandler()->ElementLocalOnBoundary(DEGREE,col);

#pragma omp for
	    for (int i=0; i<q.size(); i++)
	      {
		int iq  = q[i];
		int ile = l[i];
		
		Transformation(T,iq);
		finiteelement.ReInit(T);
	      
		GlobalToLocal(__U,u,iq);
		GlobalToLocalData(iq, __QN, __QC);
		
		integrator.BoundaryForm(*BEQ,__F,finiteelement,__U,ile,col,__QN,__QC);
		LocalToGlobal(f,__F,iq,d);
	      }
	  }
	delete BEQ;
      }
    }
    void BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const ProblemDescriptorInterface& PD, double d) const
    {
      // Do we have a boundary equation?
      if (PD.NewBoundaryEquation() == NULL)
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
      

#pragma omp parallel private(T,finiteelement,integrator,__U,__QN,__QC,__E)
      {
	auto *BEQ = PD.NewBoundaryEquation();
	const auto COLS = PD.GetBoundaryManager()->GetBoundaryEquationColors();
	BEQ->SetParameterData(QP);
      
	for(const auto col : COLS)
	  {
	    const IntVector& q = *GetDofHandler()->ElementOnBoundary(DEGREE,col);
	    const IntVector& l = *GetDofHandler()->ElementLocalOnBoundary(DEGREE,col);
#pragma omp for
	    for (int i=0; i<q.size(); i++)
	      {
		int iq  = q[i];
		int ile = l[i];
		
		Transformation(T,iq);
		finiteelement.ReInit(T);
		
		GlobalToLocal(__U,u,iq);
		GlobalToLocalData(iq, __QN, __QC);
		
		integrator.BoundaryMatrix(*BEQ,__E,finiteelement,__U,ile,col,__QN,__QC);
		LocalToGlobal(A,__E,iq,d);
	      }
	  }
	delete BEQ;
      }
    }
    double ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const 
    {
      LocalParameterData QP;
      GlobalToGlobalData(QP);

      nmatrix<double> T;

      FINITEELEMENT finiteelement;
      INTEGRATOR integrator;
      integrator.BasicInit();
      LocalVector __U;
      LocalData __QN, __QC;
      
      double j=0.;
      for(const auto col : Colors)
	{
	  const IntVector& q = *GetDofHandler()->ElementOnBoundary(DEGREE,col);
	  const IntVector& l = *GetDofHandler()->ElementLocalOnBoundary(DEGREE,col);
	  for (int i=0; i<q.size(); i++)
	    {
	      int iq  = q[i];
	      int ile = l[i];
	      
	      Transformation(T,iq);
	      finiteelement.ReInit(T);
	      
	      GlobalToLocal(__U,u,iq);
	      GlobalToLocalData(iq,__QN,__QC);
	      
	      j += integrator.ComputeBoundaryFunctional(BF,finiteelement,ile,col,__U,__QN,__QC);
	    }
	}
      
      return j;
    }
	
    ////////////////////////////////////////////////// Functionals
    double ComputePointValue(const GlobalVector& u, const Vertex2d& p0,int comp) const 
    {
      // very simple version. Only finds nodes
      for (int n=0;n<GetDofHandler()->nnodes();++n)
	{
	  double dist = 0;
	  for (int d=0;d<DIM;++d)
	    dist += pow(GetDofHandler()->vertex(n)[d]-p0[d],2.0);
	  if (dist< sqrt(1.e-13))
	    return u(n,comp);
	}
      std::cerr << "DofHandler::ComputePointValue. Vertex " << p0 << " not found!"
		<< std::endl;
      abort();
    }
    double ComputePointValue(const GlobalVector& u, const Vertex3d& p0,int comp) const 
    {
      // very simple version. Only finds nodes
      for (int n=0;n<GetDofHandler()->nnodes();++n)
	{
	  double dist = 0;
	  for (int d=0;d<DIM;++d)
	    dist += pow(GetDofHandler()->vertex(n)[d]-p0[d],2.0);
	  if (dist< sqrt(1.e-13))
	    return u(n,comp);
	}
      std::cerr << "DofHandler::ComputePointValue. Vertex " << p0 << " not found!"
		<< std::endl;
      abort();
    }
    double ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const
    {
      LocalParameterData QP;
      GlobalToGlobalData(QP);
      FP.SetParameterData(QP);
      
      int dim = GetDofHandler()->dimension();
      assert(dim==DIM);
      std::vector<int> comps = FP.GetComps();
      int nn = comps.size();
      
      std::vector<double> up(nn,0);
      
      if (dim == 2)
	{
	  auto v2d = FP.GetPoints2d();
	  assert(nn==v2d.size());
	  
	  for(int i=0;i<nn;++i)
	    {
	      up[i] = ComputePointValue(u,v2d[i],comps[i]);
	    }
	}
      else if (dim == 3)
	{
	  auto v3d = FP.GetPoints3d();
	  assert(nn==v3d.size());
	  for(int i=0;i<nn;++i)
	    {
	      up[i] = ComputePointValue(u,v3d[i],comps[i]);
	    }
	}
      else
	{
	  std::cout << "wronng dimension: dim = " << dim << std::endl;
	  abort();
	}
      
      return FP.J(up);
    }
    double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const
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
      double j=0.;
      for(int iq=0;iq<GetDofHandler()->nelements(DEGREE);++iq)
	{
	  Transformation(T,iq);
	  finiteelement.ReInit(T);
	  
	  GlobalToLocal(__U,u,iq);
	  GlobalToLocalData(iq, __QN, __QC);
	  F.point_cell(GetDofHandler()->material(DEGREE,iq));
	  j += integrator.ComputeDomainFunctional(F,finiteelement,__U,__QN,__QC);
	}
      return j;
    }
    
    
    
    ////////////////////////////////////////////////// Pressure Filter, set
    /// average zero
    void InitFilter(nvector<double> &F) const
    {
      PressureFilter *PF = static_cast<PressureFilter *>(&F);
      assert(PF);
      assert(!PF->Active());
    }


    ////////////////////////////////////////////////// Dirichlet Data
    void StrongDirichletVector(GlobalVector &u,
                               const DirichletData &BF,
                               int col,
                               const std::vector<int> &comp,
                               double d) const
    {
      const IntVector &bv = *GetDofHandler()->VertexOnBoundary(col);

      LocalParameterData QP;
      GlobalToGlobalData(QP);
      BF.SetParameterData(QP);


      FemData QH;
      nvector<double> ff(u.ncomp(), 0.);
#pragma omp parallel for private(QH) firstprivate(ff)
      for (int i = 0; i < bv.size(); ++i)
      {
        int index = bv[i];
        QH.clear();
        auto p = GetDataContainer().GetNodeData().begin();
        for (; p != GetDataContainer().GetNodeData().end(); p++)
        {
          QH[p->first].resize(p->second->ncomp());
          for (int c = 0; c < p->second->ncomp(); c++)
          {
            QH[p->first][c].m() = p->second->operator()(index, c);
          }
        }
        const Vertex<DIM> &v = GetDofHandler()->vertex(index);
        BF(ff, v, col);
        for (int iii = 0; iii < comp.size(); iii++)
        {
          int c = comp[iii];
          u(index, c) = d * ff[c];
        }
      }
    }
    void StrongDirichletMatrix(MatrixInterface &A,
                               int col,
                               const std::vector<int> &comp) const
    {
      const IntVector &bv = *GetDofHandler()->VertexOnBoundary(col);
#pragma omp parallel for
      for (int i = 0; i < bv.size(); i++)
        A.dirichlet(bv[i], comp);
    }
    void StrongDirichletMatrixOnlyRow(MatrixInterface &A,
                                      int col,
                                      const std::vector<int> &comp) const
    {
      const IntVector &bv = *GetDofHandler()->VertexOnBoundary(col);
#pragma omp parallel for
      for (int i = 0; i < bv.size(); i++)
        A.dirichlet_only_row(bv[i], comp);
    }
    void StrongDirichletVectorZero(GlobalVector &u,
                                   int col,
                                   const std::vector<int> &comp) const
    {
      const IntVector &bv = *GetDofHandler()->VertexOnBoundary(col);
#pragma omp parallel for
      for (int i = 0; i < bv.size(); i++)
      {
        int index = bv[i];
        for (int iii = 0; iii < comp.size(); iii++)
          u(index, comp[iii]) = 0.;
      }
    }


    ////////////////////////////////////////////////// Errors
    void ComputeError(const GlobalVector &u,
                      LocalVector &err,
                      const ExactSolution *ES) const
    {
      int ncomp = u.ncomp();
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
      for (int iq = 0; iq < GetDofHandler()->nelements(DEGREE); iq++)
      {
        Transformation(T, iq);
        finiteelement.ReInit(T);
        GlobalToLocal(U, u, iq);
        GlobalToLocalData(iq, QN, QC);
        integrator.ErrorsByExactSolution(lerr, finiteelement, *ES, U, QN, QC);

        for (int c = 0; c < ncomp; c++)
        {
          err(0, c) += lerr(0, c);
          err(1, c) += lerr(1, c);
          err(2, c) = std::max(err(2, c), lerr(2, c));
        }
      }
      for (int c = 0; c < ncomp; c++)
      {
        err(0, c) = sqrt(err(0, c));
        err(1, c) = sqrt(err(1, c));
      }
    }
    
    
    void ElementColoring(int degree) 
    {
      // convert mesh to graph
      std::vector<std::vector<int> > node2patch(GetDofHandler()->nnodes());

      for (int p=0;p<GetDofHandler()->nelements(degree);++p)
		{
		  const auto &iop = GetDofHandler()->GetElement(degree,p);
		  for (auto it : iop)
		    node2patch[it].push_back(p);
		}
      std::vector<std::vector<int> > adj(GetDofHandler()->nelements(degree));
      for (int p=0;p<GetDofHandler()->nelements(degree);++p)
		{
		  std::set<int> neigh;
		  const auto &iop =GetDofHandler()->GetElement(degree,p);
		  for (auto it : iop)
		    neigh.insert(node2patch[it].begin(),node2patch[it].end());
		  for (auto it : neigh)
		    if (it!=p)
		      adj[p].push_back(it);
		}
      /////////
      std::vector<int> patch2color(GetDofHandler()->nelements(degree),0);
      bool done = true;
      int color = 0;
      //faerben aller fluid patches
      do
		{
		  done = true;
		  for (int p=0;p<GetDofHandler()->nelements(degree);++p)
		    if ((patch2color[p]==color) ) // not given a color and fluid
		      {
		    done = false;
		    // block neighbors
		    for (auto n : adj[p])
		      {
		        if ((patch2color[n] == color) ) // neighbor not set, block
		          patch2color[n] = color+1;
		      }
		      }
		  ++color;
      } while (!done);
      
     
	  
      Coloring.clear();
      Coloring.resize(color-1);
      for (int p=0;p<GetDofHandler()->nelements(degree);++p)
		{
		  int col = patch2color[p];
		  assert(col<Coloring.size());
		  Coloring[col].push_back(p);
		}
      std::cout << "Coloring:" << std::endl;
    
      int colnumb=0;
      for (auto it : Coloring)
   		{
   			 std::cout<<" col " <<colnumb<< " with "<<it.size()<<"elements"<<std::endl;
   			 colnumb++;
   		}
      std::cout << std::endl;

    }
    int NumberofColors() const
	{
		return Coloring.size();
	}	

	const std::vector<int>&   elementswithcolor(int col) const
	{
	  return Coloring[col];
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

#define CGDiscQ12d                                                             \
  CGDisc<2,                                                                    \
         1,                                                                    \
         FiniteElement<2, 1, Transformation2d<BaseQ12d>, BaseQ12d>,            \
         ElementIntegratorQ12d>
#define CGDiscQ22d                                                             \
  CGDisc<2,                                                                    \
         2,                                                                    \
         FiniteElement<2, 1, Transformation2d<BaseQ22d>, BaseQ22d>,            \
         ElementIntegratorQ22d>
#define CGDiscQ13d                                                             \
  CGDisc<3,                                                                    \
         1,                                                                    \
         FiniteElement<3, 2, Transformation3d<BaseQ13d>, BaseQ13d>,            \
         ElementIntegratorQ13d>
#define CGDiscQ23d                                                             \
  CGDisc<3,                                                                    \
         2,                                                                    \
         FiniteElement<3, 2, Transformation3d<BaseQ23d>, BaseQ23d>,            \
         ElementIntegratorQ23d>


  ////// LPS
#define CGDiscQ12dLps                                                          \
  CGDisc<2,                                                                    \
         2,                                                                    \
         FiniteElement<2, 1, Transformation2d<BaseQ12dPatch>, BaseQ12dPatch>,            \
         ElementLpsIntegratorQ12d>
  
#define CGDiscQ22dLps                                                          \
  CGDisc<2,                                                                    \
         2,                                                                    \
         FiniteElement<2, 1, Transformation2d<BaseQ22d>, BaseQ22d>,            \
         ElementLpsIntegratorQ22d>

#define CGDiscQ13dLps                                                          \
  CGDisc<3,                                                                    \
         2,                                                                    \
         FiniteElement<3, 2, Transformation3d<BaseQ13dPatch>, BaseQ13dPatch>,  \
         ElementLpsIntegratorQ13d>
#define CGDiscQ23dLps                                                          \
  CGDisc<3,                                                                    \
         2,                                                                    \
         FiniteElement<3, 2, Transformation3d<BaseQ23d>, BaseQ23d>,            \
         ElementLpsIntegratorQ23d>


} // namespace Gascoigne


/*----------------------------   cgdisc.h     ---------------------------*/
/* end of #ifndef __cgdisc_H */
#endif
/*----------------------------   cgdisc.h     ---------------------------*/