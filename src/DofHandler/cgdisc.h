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

#include  "discretizationinterface.h"
#include  "sparsestructure.h"
#include  "mginterpolatornested.h"
#include  "pressurefilter.h"

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
  template<int DIM, int DEGREE, class FINITEELEMENT, class INTEGRATOR>
  class CGDisc : public DiscretizationInterface
  {
    private:

    const DofHandler<DIM>*  dofhandler;
    DataContainer           datacontainer;

    protected:


    public:
    
    CGDisc(){}
    ~CGDisc(){}

    const DofHandler<DIM>* GetDofHandler() const { return dofhandler; }
    

    const DataContainer& GetDataContainer() const { return datacontainer; }
    void SetDataContainer(const DataContainer& q) { datacontainer = q; }

    

      //
      //// Functions called from the Solver
      //
    std::string GetName() const { return "CG Discretization";}

    void AddNodeVector(const std::string& name, const GlobalVector* q) const {assert(0);}
    void DeleteNodeVector(const std::string& name) const {assert(0);}

    void AddCellVector(const std::string& name, const GlobalVector* q) const {assert(0);}
    void DeleteCellVector(const std::string& name) const {assert(0);}
    
    void AddParameterVector(const std::string& name, const GlobalParameterVector* q) const {assert(0);}
    void DeleteParameterVector(const std::string& name) const {assert(0);}

    void BasicInit(const ParamFile* pf)
    {
      // HANGING NODES
    }
    void ReInit   (const GascoigneMesh* M)
    {
      dofhandler = dynamic_cast<const DofHandler<DIM>*> (M);
      assert(dofhandler);
      // HANGING NODES
    }

    int n() const { return GetDofHandler()->nnodes(); }
    int nc() const { return GetDofHandler()->nelements(DEGREE); }
    int n_withouthanging()const {
      return n();
      // HANGING NODES
    }
    

    //////////////////////////////////////////////////
    void Transformation(nmatrix<double>& T, int iq) const
    {
      assert (GetDofHandler()->dimension() == DIM);
      int ne = GetDofHandler()->nodes_per_element(DEGREE);
      
      IntVector indices = GetDofHandler()->GetElement(DEGREE,iq);
      assert(ne==indices.size());
      
      T.memory(DIM,ne);
      if(DIM==2)
	{
	  for(int ii=0;ii<ne;ii++)
	    {
	      Vertex2d v = GetDofHandler()->vertex2d(indices[ii]);
	      T(0,ii) = v.x();               
	      T(1,ii) = v.y();
	    }
	}
      else if(DIM==3)
	{
	  for(int ii=0;ii<ne;ii++)
	    {
	      Vertex3d v = GetDofHandler()->vertex3d(indices[ii]);
	      T(0,ii) = v.x();               
	      T(1,ii) = v.y();
	      T(2,ii) = v.z();
	    }
	}
    }
    void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
    {
      MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
      assert(IP);
      IP->BasicInit(MT);
    }

    
    void Structure(SparseStructureInterface* SI) const
    {
      SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
      assert(S);
      
      S->build_begin(n());
      for(int iq=0;iq<GetDofHandler()->nelements(DEGREE);iq++)
	{
	  IntVector indices = GetDofHandler()->GetElement(DEGREE,iq);
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
      
      for(auto p : gpd)
	QP.insert(make_pair(p.first,*p.second));
    }
    void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
    {
      IntVector indices = GetDofHandler()->GetElement(DEGREE,iq);
      
      // HANGING NODES
      //HN->CondenseHanging(E,indices);
      IntVector::const_iterator  start = indices.begin();
      IntVector::const_iterator  stop  = indices.end();
      A.entry(start,stop,E,s);
    }
    void LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const
    {
      IntVector indices = GetDofHandler()->GetElement(DEGREE,iq);
      for(int ii=0; ii<indices.size(); ii++) 
	{
	  int i = indices[ii];
	  f.add_node(i,s,ii,F);
	}
    }
    void GlobalToLocal(LocalVector& U, const GlobalVector& u, int iq) const
    {
      IntVector indices = GetDofHandler()->GetElement(DEGREE,iq);
      U.ReInit(u.ncomp(),indices.size());
      for(int ii=0; ii<indices.size(); ii++) 
	{
	  int i = indices[ii];
	  U.equ_node(ii,i,u);
	}
    }
    void GlobalToLocalCell(LocalVector& U, const GlobalVector& u, int iq) const
    {
      U.ReInit(u.ncomp(),1);
      for(int c=0;c<u.ncomp();++c)
	U(0,c) = u(iq,c);
    }
    void GlobalToLocalData(int iq, LocalData QN, LocalData QC) const
    {
      const GlobalData& gnd = GetDataContainer().GetNodeData();
      QN.clear();
      GlobalData::const_iterator p=gnd.begin();
      for(; p!=gnd.end(); p++)
	{
	  GlobalToLocal(QN[p->first],*p->second,iq);
	}
      
      const GlobalData& gcd = GetDataContainer().GetCellData();
      QC.clear();
      GlobalData::const_iterator q=gcd.begin();
      for(; q!=gcd.end(); q++)
	{
	  GlobalToLocalCell(QC[q->first],*q->second,iq);
	}
    }
    
    


    ////////////////////////////////////////////////// integration ueber Zellen
    void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
    {
      LocalParameterData QP;
      GlobalToGlobalData(QP);
      EQ.SetParameterData(QP);

      nmatrix<double> T;
      FINITEELEMENT   finiteelement;
      INTEGRATOR      integrator;
      integrator.BasicInit();
      LocalVector     __U,__F;
      LocalData       __QN, __QC;
#pragma omp parallel for private(T,finiteelement,integrator,__U,__F,__QN, __QC) 
      for(int iq=0;iq<GetDofHandler()->nelements(DEGREE);++iq)
	{
	  Transformation(T,iq);
	  finiteelement.ReInit(T);
	  
	  GlobalToLocal(__U,u,iq);
	  GlobalToLocalData(iq,__QN,__QC);
	  EQ.point_cell(GetDofHandler()->material(DEGREE,iq));
	  
	  integrator.Form(EQ,__F,finiteelement,__U,__QN,__QC);
	  LocalToGlobal(f,__F,iq,d);
	}
      
    }
    void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
    {
      LocalParameterData QP;
      GlobalToGlobalData(QP);
      RHS.SetParameterData(QP);
      

      nmatrix<double> T;
      FINITEELEMENT   finiteelement;
      INTEGRATOR      integrator;
      integrator.BasicInit();
      LocalVector     __F;
      LocalData       __QN, __QC;
      
      for(int iq=0;iq<GetDofHandler()->nelements(DEGREE);++iq)
	{
	  Transformation(T,iq);
	  finiteelement.ReInit(T);
	  
	  GlobalToLocalData(iq, __QN, __QC);
	  RHS.point_cell(GetDofHandler()->material(DEGREE,iq));
	  integrator.Rhs(RHS,__F,finiteelement,__QN,__QC);
	  LocalToGlobal(f,__F,iq,s);
	}
    }
    void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
    {
      LocalParameterData QP;
      GlobalToGlobalData(QP);
      EQ.SetParameterData(QP);
      
      nmatrix<double> T;
      FINITEELEMENT   finiteelement;
      INTEGRATOR      integrator;
      integrator.BasicInit();
      LocalVector     __F, __U;
      LocalData       __QN, __QC;
      EntryMatrix     __E;
      
      for(int iq=0;iq<GetDofHandler()->nelements(DEGREE);++iq)
	{
	  Transformation(T,iq);
	  finiteelement.ReInit(T);

	  GlobalToLocal(__U,u,iq);
	  GlobalToLocalData(iq, __QN, __QC);
	  
	  EQ.point_cell(GetDofHandler()->material(DEGREE, iq));
	  //EQ.cell(GetDofHandler(),iq,__U,__QN);
	  integrator.Matrix(EQ,__E,finiteelement,__U,__QN,__QC);
	  LocalToGlobal(A,__E,iq,d);
	}
      // HANGING NODES
      //    HN->MatrixDiag(u.ncomp(),A);
    }

    
    ////////////////////////////////////////////////// Pressure Filter, set average zero
    void InitFilter(nvector<double>& F) const
    {
      PressureFilter* PF = static_cast<PressureFilter*>(&F);
      assert(PF);
      assert (!PF->Active());
    }




    ////////////////////////////////////////////////// Dirichlet Data
    void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp, double d) const
    {
      nvector<double> ff(u.ncomp(),0.);
      const IntVector& bv = *GetDofHandler()->VertexOnBoundary(col);
      
      LocalParameterData  QP;
      GlobalToGlobalData(QP);
      BF.SetParameterData(QP);


      FemData QH;
      for(auto index : bv)
	{
	  QH.clear();
	  auto p = GetDataContainer().GetNodeData().begin();
	  for(; p!=GetDataContainer().GetNodeData().end(); p++)
	    {
	      QH[p->first].resize(p->second->ncomp());
	      for(int c=0; c<p->second->ncomp(); c++)
		{
		  QH[p->first][c].m() = p->second->operator()(index,c);
		}
	    }
	  const Vertex<DIM>& v = GetDofHandler()->vertex(index);
	  BF(ff,v,col);
	  for(int iii=0;iii<comp.size();iii++)
	    {
	      int c = comp[iii];
	      u(index,c) = d * ff[c];
	    }
	}
    }
    void StrongDirichletMatrix(MatrixInterface& A, int col, const std::vector<int>& comp) const
    {
      const IntVector& bv = *GetDofHandler()->VertexOnBoundary(col);
      for(int i=0;i<bv.size();i++)
	A.dirichlet(bv[i], comp);
    }
    void StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const std::vector<int>& comp) const
    {
      const IntVector& bv = *GetDofHandler()->VertexOnBoundary(col);
      for(int i=0;i<bv.size();i++)
	A.dirichlet_only_row(bv[i], comp);
    }   
    void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const
    {
      const IntVector& bv = *GetDofHandler()->VertexOnBoundary(col);
      for(int i=0;i<bv.size();i++)
	{
	  int index = bv[i];
	  for(int iii=0;iii<comp.size();iii++)
	    u( index,comp[iii] ) = 0.;
	}  
    }




    ////////////////////////////////////////////////// Errors
    void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const
    {
      int ncomp = u.ncomp();
      err.ncomp() = ncomp;
      err.reservesize(3);
      err = 0.;
      
      GlobalVector lerr(ncomp,3); 
      lerr.zero();
      


      LocalParameterData  QP;
      GlobalToGlobalData(QP);
      ES->SetParameterData(QP);
      

      nmatrix<double> T;
      LocalVector     U;
      LocalData       QN,QC;
      FINITEELEMENT finiteelement;
      INTEGRATOR    integrator;
      integrator.BasicInit();
      for(int iq=0; iq<GetDofHandler()->nelements(DEGREE); iq++)
	{
	  Transformation(T,iq);
	  finiteelement.ReInit(T);
	  GlobalToLocal(U,u,iq);
	  GlobalToLocalData(iq,QN,QC);
	  integrator.ErrorsByExactSolution(lerr,finiteelement,*ES,U,QN,QC);
	  
	  for(int c=0;c<ncomp;c++)  
	    {
	      err(0,c) += lerr(0,c);
	      err(1,c) += lerr(1,c);
	      err(2,c) = std::max(err(2,c),lerr(2,c));
	    }
	}
      for(int c=0;c<ncomp;c++)  
	{
	  err(0,c) = sqrt(err(0,c));
	  err(1,c) = sqrt(err(1,c));
	}
    }
  };
}


/*----------------------------   cgdisc.h     ---------------------------*/
/* end of #ifndef __cgdisc_H */
#endif
/*----------------------------   cgdisc.h     ---------------------------*/
