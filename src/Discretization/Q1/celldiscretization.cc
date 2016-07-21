/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2011 by the Gascoigne 3D authors
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


#include  "celldiscretization.h"
#include  "sparsestructure.h"
#include  "pressurefilter.h"
#include  "gascoignemesh.h"
#include  "columndiagstencil.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
#include  "simplematrix.h"
#include  "stopwatch.h"
#include  "hnstructureq12d.h"
#include  "hnstructureq13d.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "gascoignemeshtransfer.h"
#include "galerkinintegrator.h"
#include <omp.h>

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{

  ///////////// Initialization
  
  template<int DIM>
  void Q1<DIM>::BasicInit(const ParamFile* pf)
  {
    assert(HN==NULL);
    HN = NewHNStructure();
    assert(HN);
    
    if(!GetIntegratorPointer())
      GetIntegratorPointer() =  new GalerkinIntegrator<DIM>;
    assert(GetIntegrator());
    
    GetIntegratorPointer()->BasicInit();
    
    if(!GetFemPointer())
      {
	//      typedef Transformation2d<BaseQ1<2> >           TransQ1;
	typedef Gascoigne::Transformation<DIM, BaseQ1<DIM> >   TransQ1;
	typedef FiniteElement<DIM,DIM-1,TransQ1,BaseQ1<DIM> >  FiniteElement;
	Q1<DIM>::GetFemPointer()   =  new FiniteElement;
      }
    assert(GetFem());
  }
  


  template<int DIM>
  void Q1<DIM>::InitColoring() 
  {
    RealTimeStopWatch rt;
    rt.start();
    
    //// make graph for coloring
    // first node to cell
    vector<set<int> > node2cell(GetMesh()->nnodes());
    for (int c=0;c<GetMesh()->ncells();++c)
      {
	IntVector indices = GetLocalIndices(c);
	for (int j=0;j<indices.size();++j) 
	  node2cell[indices[j]].insert(c);
      }
    
    vector<int> neighbors;
    vector<int> start;
    int index =0;
    for (int c=0;c<GetMesh()->ncells();++c)
      {
	start.push_back(index);
	IntVector indices = GetLocalIndices(c);
	set<int> n; // neighbors of cell c
	for (int j=0;j<indices.size();++j) 
	  for (set<int>::const_iterator it =   node2cell[indices[j]].begin();
	       it!=node2cell[indices[j]].end();++it)
	    n.insert(*it);
	for (set<int>::const_iterator it = n.begin();it!=n.end();++it)
	  {
	    if (*it!=c) 
	      {
		neighbors.push_back(*it);
		++index;
	      }
	  }
      }
    start.push_back(index);

    
    // partition graph
    assert(GetMesh()->ncells()+1==start.size());
    vector<int> cell2color(GetMesh()->ncells(),-1);
    // ganz primitiv, immer kleinste freie Nummer suchen.
    for (int c=0;c<GetMesh()->ncells();++c)
      {
	set<int> nc; // farben der nachbarn
	for (int ni=start[c];ni<start[c+1];++ni)
	  nc.insert(cell2color[neighbors[ni]]);
	int col = 0;
	while (nc.find(col)!=nc.end()) ++col;
	cell2color[c]=col;
      }

    _col_graph_ncol=0;
    for (int c=0;c<cell2color.size();++c)
      _col_graph_ncol = std::max(_col_graph_ncol,cell2color[c]);
    _col_graph_ncol++;
    _col_graph.resize(_col_graph_ncol);
    for (int c=0;c<cell2color.size();++c)
      _col_graph[cell2color[c]].push_back(c);
    
    // Statistics
    /*
      vector<int> histo(20);
      for (int c=0;c<GetMesh()->ncells();++c)
      {
      assert(cell2color[c]<20);
      histo[cell2color[c]]++;
      }
      for (int i=0;i<20;++i)
      cout << i << " " << histo[i] << endl;
    */
    rt.stop();
    /*
      std::cout << "Coloring Graph " << neighbors.size() << " " << start.size() << std::endl;
      std::cout << "Coloring time " << rt.read() << std::endl;
    */
  }
  




  template<int DIM>
  void Q1<DIM>::Structure(SparseStructureInterface* SI) const
  {
    SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
    assert(S);
    
    S->build_begin(n());
    for(int iq=0;iq<GetMesh()->ncells();iq++)
      {
	IntVector indices = GetLocalIndices(iq);
	HN->CondenseHanging(indices);
	S->build_add(indices.begin(), indices.end());
      }
    HN->SparseStructureDiag(S);
    S->build_end();  
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::Transformation(FemInterface::Matrix& T, int iq) const
  {
    int dim = GetMesh()->dimension();
    int ne = GetMesh()->nodes_per_cell(iq);
	
    IntVector indices = GetMesh()->IndicesOfCell(iq);
    assert(ne==indices.size());
	
    T.memory(dim,ne);
    if(dim==2)
      {
	for(int ii=0;ii<ne;ii++)
	  {
	    Vertex2d v = GetMesh()->vertex2d(indices[ii]);
	    T(0,ii) = v.x();               
	    T(1,ii) = v.y();
	  }
      }
    else if(dim==3)
      {
	for(int ii=0;ii<ne;ii++)
	  {
	    Vertex3d v = GetMesh()->vertex3d(indices[ii]);
	    T(0,ii) = v.x();               
	    T(1,ii) = v.y();
	    T(2,ii) = v.z();
	  }
      }
  }

  /* ----------------------------------------- */
  /* ----------------------------------------- */
  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
  {
    nmatrix<double> T;
  
    GlobalToGlobalData();
    EQ.SetParameterData(__QP);

    assert(_col_graph_ncol>0);

    LocalVector _U=__U,_F=__F;
    
    map<const string,LocalVector> _QN;

    typename FemType<DIM>::FFF FEM;
    
    vector<GalerkinIntegrator<DIM> > INT(8);
    for (int i=0;i<INT.size();++i)
      INT[i].BasicInit();
    
    omp_set_nested(1);
    
    for (int col=0;col<_col_graph_ncol;++col)
      {
	//#pragma omp parallel for num_threads(8) private(T,_U,_F,FEM,_QN)
	for(int ii=0;ii<_col_graph[col].size();++ii)
	  {
	    int iq = _col_graph[col][ii];
	    Transformation(T,iq);
	    //GetFem()->ReInit(T);
	    FEM.ReInit(T);
	    
	    //	    GlobalToLocal(_U,u,iq);
	    GlobalToLocalSingle(_U,u,iq);
	    _QN.clear();
	    for (GlobalData::const_iterator p=GetDataContainer().GetNodeData().begin();
		 p!=GetDataContainer().GetNodeData().end(); p++)
	      GlobalToLocalSingle(_QN[p->first],*p->second,iq);
	    //
	    
	    assert(__QC.size()==0);

	    assert(omp_get_thread_num()<INT.size());
	    //INT[omp_get_thread_num()].
	    GetIntegrator()->Form(EQ,_F,FEM,_U,_QN,__QC);
	    
	    BasicDiscretization::LocalToGlobal(f,_F,iq,d);
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
  {
    nmatrix<double> T;
  
    GlobalToGlobalData();
    EQ.SetParameterData(__QP);
  
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GlobalToLocal(__U,u,iq);
	//EQ.cell(GetMesh(),iq,__U,__QN);
	GetIntegrator()->AdjointForm(EQ,__F,*GetFem(),__U,__QN,__QC);
	BasicDiscretization::LocalToGlobal(f,__F,iq,d);
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
  {
    nmatrix<double> T;
  
    GlobalToGlobalData();
    BE.SetParameterData(__QP);
  
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];

	    Transformation(T,iq);
	    GetFem()->ReInit(T);

	    GlobalToLocal(__U,u,iq);

	    GetIntegrator()->BoundaryForm(BE,__F,*GetFem(),__U,ile,col,__QN,__QC);
	    BasicDiscretization::LocalToGlobal(f,__F,iq,d);
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::Matrix(MatrixInterface& A, const GlobalVector& u, const ProblemDescriptorInterface* PD, double d) const
  {
    const Equation& EQ = *(PD->GetEquation());
  
    nmatrix<double> T;
  
    GlobalToGlobalData();
    EQ.SetParameterData(__QP);
  
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GlobalToLocal(__U,u,iq);
	//EQ.cell(GetMesh(),iq,__U,__QN);
	GetIntegrator()->Matrix(EQ,__E,*GetFem(),__U,__QN,__QC);
	LocalToGlobal(A,__E,iq,d);
      }

    HN->MatrixDiag(u.ncomp(),A);
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
  {
    nmatrix<double> T;
  
    GlobalToGlobalData();
    BE.SetParameterData(__QP);
  
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];
          
	    Transformation(T,iq);
	    GetFem()->ReInit(T);

	    GlobalToLocal(__U,u,iq);
	    GetIntegrator()->BoundaryMatrix(BE,__E,*GetFem(),__U,ile,col,__QN,__QC);
	    LocalToGlobal(A,__E,iq,d);
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::MassMatrix(MatrixInterface& A) const
  {
    nmatrix<double> T;
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);
	GetIntegrator()->MassMatrix(__E,*GetFem());
	LocalToGlobal(A,__E,iq,1.);
      }
    HN->MatrixDiag(1,A);  
  }

  /* ------------------------------------------ */

  template<int DIM>
  void Q1<DIM>::BoundaryMassMatrix(MatrixInterface& A, const IntSet& Colors) const
  {
    A.zero();
    nmatrix<double> T;
  
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];
          
	    Transformation(T,iq);
	    GetFem()->ReInit(T);

	    GetIntegrator()->BoundaryMassMatrix(__E,*GetFem(),ile);
	    LocalToGlobal(A,__E,iq,1.);
	  }
      }
    //Diagonaleintraege auf 1 setzen, wenn Eintrag noch null, damit A invertierbar ist.
    const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil*>(A.GetStencil());
    assert(ST);
    SimpleMatrix* SM = dynamic_cast<SimpleMatrix*>(&A);
    assert(SM);
    int n = ST->n();
    int pos;
    for(int i = 0; i < n; i++)
      {
	pos = ST->diag(i);
	if(SM->GetValue(pos) == 0)
	  {
	    SM->GetValue(pos) = 1;
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Gascoigne::Q1<DIM>::MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const
  {
    nmatrix<double> T;
 
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GlobalToLocal(__U,u,iq);
	GetIntegrator()->MassForm(TP,__F,*GetFem(),__U);
	BasicDiscretization::LocalToGlobal(f,__F,iq,s);
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const
  {
    int ncomp = u.ncomp();
    err.ncomp() = ncomp;
    err.reservesize(3);
    err = 0.;

    GlobalVector lerr(ncomp,3); 
    lerr.zero();

    nmatrix<double> T;
  
    GlobalToGlobalData();
    ES->SetParameterData(__QP);
  
    for(int iq=0; iq<GetMesh()->ncells(); iq++)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);
	GlobalToLocal(__U,u,iq);
	GetIntegrator()->ErrorsByExactSolution(lerr,*GetFem(),*ES,__U,__QN,__QC);

	for(int c=0;c<ncomp;c++)  
	  {
	    err(0,c) += lerr(0,c);
	    err(1,c) += lerr(1,c);
	    err(2,c) = Gascoigne::max(err(2,c),lerr(2,c));
	  }
      }
    for(int c=0;c<ncomp;c++)  
      {
	err(0,c) = sqrt(err(0,c));
	err(1,c) = sqrt(err(1,c));
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::AssembleError(GlobalVector& eta, const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const
  {
    int ncomp = u.ncomp();
    err.ncomp() = ncomp;
    err.reservesize(3);
    err = 0.;

    eta.ncomp() = 3;
    eta.reservesize(GetMesh()->ncells());

    GlobalVector lerr(ncomp,3); 
    lerr.zero();

    nmatrix<double> T;
  
    GlobalToGlobalData();
    ES->SetParameterData(__QP);
  
    for(int iq=0; iq<GetMesh()->ncells(); iq++)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);
	GlobalToLocal(__U,u,iq);
	GetIntegrator()->ErrorsByExactSolution(lerr,*GetFem(),*ES,__U,__QN,__QC);

	for(int c=0;c<ncomp;c++)  
	  {
	    err(0,c) += lerr(0,c);
	    err(1,c) += lerr(1,c);
	    err(2,c) = Gascoigne::max(err(2,c),lerr(2,c));

	    eta(iq,0) += lerr(0,c);
	    eta(iq,1) += lerr(1,c);
	    eta(iq,2) += lerr(2,c);
	  }
	eta(iq,0) = sqrt(eta(iq,0));
	eta(iq,1) = sqrt(eta(iq,1));
	eta(iq,2) = sqrt(eta(iq,2));
      }
    for(int c=0;c<ncomp;c++)  
      {
	err(0,c) = sqrt(err(0,c));
	err(1,c) = sqrt(err(1,c));
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
  {
    nmatrix<double> T;
  
    GlobalToGlobalData();
    RHS.SetParameterData(__QP);
  
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GlobalToLocalData(iq);
	GetIntegrator()->Rhs(RHS,__F,*GetFem(),__QN,__QC);
	BasicDiscretization::LocalToGlobal(f,__F,iq,s);
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const
  {
    nmatrix<double> T;
  
    GlobalToGlobalData();
    NRHS.SetParameterData(__QP);
  
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];

	    Transformation(T,iq);
	    GetFem()->ReInit(T);

	    GlobalToLocalData(iq);
	    GetIntegrator()->BoundaryRhs(NRHS,__F,*GetFem(),ile,col,__QN,__QC);
	    BasicDiscretization::LocalToGlobal(f,__F,iq,s);
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::InitFilter(DoubleVector& F) const
  {
    PressureFilter* PF = static_cast<PressureFilter*>(&F);
    assert(PF);

    if (!PF->Active()) return;

    PF->ReInit(GetMesh()->nnodes(),GetMesh()->nhanging());
    nmatrix<double> T;
    for(int iq=0; iq<GetMesh()->ncells(); ++iq)
      {
	int nv = GetMesh()->nodes_per_cell(iq);
	EntryMatrix  E(nv,1);

	Transformation(T,iq);
	GetFem()->ReInit(T);

	double cellsize = GetIntegrator()->MassMatrix(E,*GetFem());
	PF->AddDomainPiece(cellsize);

	IntVector ind = GetMesh()->IndicesOfCell(iq);
	HN->CondenseHanging(E,ind);

	for(int i=0;i<ind.size();i++)
	  {
	    for(int j=0;j<ind.size();j++)
	      {
		F[ind[j]] += E(i,j,0,0);
	      }
	  }
      }
  }

  /* ----------------------------------------- */

  template<>
  void Q1<2>::DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const
  {
    int dim = GetMesh()->dimension();
    assert(dim==2);
    vector<int> comps = DRHS.GetComps();
    int nn = comps.size();

    vector<double> up(nn,0);
 
    vector<Vertex2d> v2d = DRHS.GetPoints2d();
    assert(nn==v2d.size());
    
    for(int i=0;i<nn;++i)
      {
	DiracRhsPoint(f,DRHS,v2d[i],i,s);
      }
  }
  template<>
  void Q1<3>::DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const
  {
    int dim = GetMesh()->dimension();
    assert(dim==3);
    vector<int> comps = DRHS.GetComps();
    int nn = comps.size();

    vector<double> up(nn,0);
 
    vector<Vertex3d> v3d = DRHS.GetPoints3d();
    assert(nn==v3d.size());
    for(int i=0;i<nn;++i)
      {
	DiracRhsPoint(f,DRHS,v3d[i],i,s);
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex<DIM>& p0,int i,double s) const
  {
    __F.ReInit(f.ncomp(),GetFem()->n());

    Vertex<DIM> Tranfo_p0;
    
    int iq = GetCellNumber(p0,Tranfo_p0);
    if (iq==-1)
      {
	cerr << "Q1<DIM>::DiracRhsPoint point not found\n";
	abort();
      }
  
    nmatrix<double> T;
    Transformation(T,iq);
    GetFem()->ReInit(T);
  
    GlobalToLocalData(iq);
    GlobalToGlobalData();
    DRHS.SetParameterData(__QP);

    GetIntegrator()->DiracRhsPoint(__F,*GetFem(),Tranfo_p0,DRHS,i,__QN,__QC);
    BasicDiscretization::LocalToGlobal(f,__F,iq,s);
  }

  /* ----------------------------------------- */


  template<int DIM>
  double Q1<DIM>::ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const 
  {
    GlobalToGlobalData();
    BF.SetParameterData(__QP);
  
    nmatrix<double> T;
    double j=0.;
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];

	    Transformation(T,iq);
	    GetFem()->ReInit(T);

	    GlobalToLocal(__U,u,iq);
	    j += GetIntegrator()->ComputeBoundaryFunctional(BF,*GetFem(),ile,col,__U,__QN,__QC);
	  }
      }

    return j;
  }

  /* ----------------------------------------- */

  template<int DIM>
  double Q1<DIM>::ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const
  {
    GlobalToGlobalData();
    F.SetParameterData(__QP);
  
    nmatrix<double> T;
    double j=0.;
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GlobalToLocal(__U,u,iq);
	j += GetIntegrator()->ComputeDomainFunctional(F,*GetFem(),__U,__QN,__QC);
      }
    return j;
  }

  /* ----------------------------------------- */

  template<int DIM>
  double Q1<DIM>::ComputeErrorDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const
  {
    GlobalToGlobalData();
    F.SetParameterData(__QP);

    nmatrix<double> T;
    double j=0.;
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GlobalToLocal(__U,u,iq);
	j += GetIntegrator()->ComputeErrorDomainFunctional(F,*GetFem(),__U,__QN,__QC);
      }
    return j;
  }

  /* ----------------------------------------- */

  template<>
  double Q1<2>::ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const
  {
    FP.SetParameterData(__QP);

    int dim = GetMesh()->dimension(); assert(dim==2);
    vector<int> comps = FP.GetComps();
    int nn = comps.size();

    vector<double> up(nn,0);
 
    vector<Vertex2d> v2d = FP.GetPoints2d();
    assert(nn==v2d.size());
    
    for(int i=0;i<nn;++i)
      up[i] = ComputePointValue(u,v2d[i],comps[i]);

    return FP.J(up);
  }
  template<>
  double Q1<3>::ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const
  {
    FP.SetParameterData(__QP);

    int dim = GetMesh()->dimension(); assert(dim==3);
    vector<int> comps = FP.GetComps();
    int nn = comps.size();

    vector<double> up(nn,0);

    vector<Vertex3d> v3d = FP.GetPoints3d();
    assert(nn==v3d.size());
    for(int i=0;i<nn;++i)
      up[i] = ComputePointValue(u,v3d[i],comps[i]);

    return FP.J(up);
  }
  /* ----------------------------------------- */

  template<int DIM>
  double Q1<DIM>::ComputePointValue(const GlobalVector& u, const Vertex<DIM>& p0,int comp) const
  {
    Vertex<DIM> Tranfo_p0;
    
    int iq = GetCellNumber(p0,Tranfo_p0);
    if (iq==-1)
      {
	cerr << "Q1<DIM>::ComputePointValue point not found\n";
	abort();
      }

    nmatrix<double> T;
    Transformation(T,iq);
    GetFem()->ReInit(T);

    GlobalToLocal(__U,u,iq);

    return GetIntegrator()->ComputePointValue(*GetFem(),Tranfo_p0,__U,comp);
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::EvaluateCellRightHandSide(GlobalVector& f, const DomainRightHandSide& CF, double d) const
  {
    nmatrix<double> T;

    GlobalToGlobalData();
    CF.SetParameterData(__QP);

    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GlobalToLocalData(iq);
	GetIntegrator()->EvaluateCellRightHandSide(__F,CF,*GetFem(),__QN,__QC);

	f.add_node(iq,d,0,__F);
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::EvaluateBoundaryCellRightHandSide(GlobalVector& f,const IntSet& Colors, const BoundaryRightHandSide& CF, double d) const
  {
    nmatrix<double> T;
  
    GlobalToGlobalData();
    CF.SetParameterData(__QP);
  
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];

	    Transformation(T,iq);
	    GetFem()->ReInit(T);

	    GlobalToLocalData(iq);
	    GetIntegrator()->EvaluateBoundaryCellRightHandSide(__F,CF,*GetFem(),ile,col,__QN,__QC);

	    f.add_node(iq,d,0,__F);
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::EvaluateParameterRightHandSide(GlobalVector& f, const DomainRightHandSide& CF, double d) const
  {
    nmatrix<double> T;

    GlobalToGlobalData();
    CF.SetParameterData(__QP);

    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GlobalToLocalData(iq);
	GetIntegrator()->EvaluateCellRightHandSide(__F,CF,*GetFem(),__QN,__QC);

	f.add(d,__F);
      }
  }

  /* ----------------------------------------- */ 
  
  template<int DIM>
  void Gascoigne::Q1<DIM>::EvaluateBoundaryParameterRightHandSide(GlobalVector& f,const IntSet& Colors,
								  const BoundaryRightHandSide& CF, double d) const
  {
    nmatrix<double> T;
  
    GlobalToGlobalData();
    CF.SetParameterData(__QP);
  
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];

	    Transformation(T,iq);
	    GetFem()->ReInit(T);

	    GlobalToLocalData(iq);
	    GetIntegrator()->EvaluateBoundaryCellRightHandSide(__F,CF,*GetFem(),ile,col,__QN,__QC);

	    f.add(d,__F);
	  }
      }
  }

  /* ----------------------------------------- */ 

  template<int DIM>  
  void Q1<DIM>::InterpolateDomainFunction(GlobalVector& f, const DomainFunction& DF) const
  {
    int dim = GetMesh()->dimension();
    f.zero();

    DoubleVector gf;
    gf.resize(DF.GetNcomp());

    const GlobalData& gnd = GetDataContainer().GetNodeData();
    FemData QH;

    if(dim==2)
      {
	for(int r=0; r<GetMesh()->nnodes(); ++r)
	  {
	    QH.clear();
	    GlobalData::const_iterator p=gnd.begin();
	    for(; p!=gnd.end(); p++)
	      {
		FemFunction& UH = QH[p->first];
		const GlobalVector& U = *p->second;
		UH.resize(U.ncomp());
		for (int c=0; c<UH.size(); c++)
		  {
		    UH[c].zero();
		    UH[c].m() = U(r,c);
		  }
	      }

	    DF.SetFemData(QH);
	    Vertex2d v = GetMesh()->vertex2d(r);
	    DF.F(gf,v);
	    f.add_node(r,1.,gf);
	  }
      }
    else
      {
	for(int r=0; r<GetMesh()->nnodes(); ++r)
	  {
	    QH.clear();
	    GlobalData::const_iterator p=gnd.begin();
	    for(; p!=gnd.end(); p++)
	      {
		FemFunction& UH = QH[p->first];
		const GlobalVector& U = *p->second;
		for (int c=0; c<UH.size(); c++)
		  {
		    UH[c].zero();
		    UH[c].m() = U(r,c);
		  }
	      }

	    DF.SetFemData(QH);
	    Vertex3d v = GetMesh()->vertex3d(r);
	    DF.F(gf,v);
	    f.add_node(r,1.,gf);
	  }
      }
  }

  /* ----------------------------------------- */  

  template<int DIM>
  void Q1<DIM>::InterpolateCellDomainFunction(GlobalVector& f, const DomainFunction& DF) const
  {
    int dim = GetMesh()->dimension();
    f.zero();

    DoubleVector gf;
    gf.resize(DF.GetNcomp());

    if(dim==2)
      {
	Vertex2d v;
	for(int iq=0;iq<GetMesh()->ncells();++iq)
	  {
	    v.zero();
	    for (int in=0; in<4; in++)
	      {
		int r = GetMesh()->vertex_of_cell(iq,in);
		v += GetMesh()->vertex2d(r);
	      }
	    v *= 0.25;
	    DF.F(gf,v);
	    f.add_node(iq,1.,gf);
	  }
      }
    else
      {
	Vertex3d v;
	for(int iq=0;iq<GetMesh()->ncells();++iq)
	  {
	    v.zero();
	    for (int in=0; in<8; in++)
	      {
		int r = GetMesh()->vertex_of_cell(iq,in);
		v += GetMesh()->vertex3d(r);
	      }
	    v *= 0.125;
	    DF.F(gf,v);
	    f.add_node(iq,1.,gf);
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::Transformation_HM(FemInterface::Matrix& T, const HierarchicalMesh* HM, int iq) const
  {
    int dim = GetMesh()->dimension();
    int ne = GetMesh()->nodes_per_cell(iq);

    IntVector indices = HM->GetVertices(iq);
    assert(ne==indices.size());
    swapIndices(indices);

    T.memory(dim,ne);
    if(dim==2)
      {
	for(int ii=0;ii<ne;ii++)
	  {
	    Vertex2d v = GetMesh()->vertex2d(indices[ii]);
	    T(0,ii) = v.x();               
	    T(1,ii) = v.y();
	  }
      }
    else if(dim==3)
      {
	for(int ii=0;ii<ne;ii++)
	  {
	    Vertex3d v = GetMesh()->vertex3d(indices[ii]);
	    T(0,ii) = v.x();               
	    T(1,ii) = v.y();
	    T(2,ii) = v.z();
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::GlobalToLocal_HM(LocalVector& U, const GlobalVector& u, const HierarchicalMesh* HM, int iq) const
  {
    IntVector indices = HM->GetVertices(iq);
    swapIndices(indices);

    U.ReInit(u.ncomp(),indices.size());
    for(int ii=0; ii<indices.size(); ii++) 
      {
	int i = indices[ii];
	U.equ_node(ii,i,u);
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::swapIndices(IntVector& indices) const
  {
    assert(indices.size()>=4);

    int help = indices[2];
    indices[2] = indices[3];
    indices[3] = help;
    if (indices.size()==8)
      {
	help = indices[6];
	indices[6] = indices[7];
	indices[7] = help;
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::GetVolumes(DoubleVector& a) const
  {
    a.resize(GetMesh()->ncells());
    nmatrix<double> T;
    int dim = GetMesh()->dimension();
  
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);
	if(dim==2)
	  {
	    Vertex2d xi;
	    xi.x() = 0.5;
	    xi.y() = 0.5;
	    GetFem()->point(xi);
	  }
	else
	  {
	    Vertex3d xi;
	    xi.x() = 0.5;
	    xi.y() = 0.5;
	    xi.z() = 0.5;
	    GetFem()->point(xi);
	  }
	a[iq] = GetFem()->J();
      }
  }
  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::GetAreas(DoubleVector& a, const IntSet& Colors) const
  {
    a.resize(GetMesh()->ncells(),1.);
    nmatrix<double> T;
    int dim = GetMesh()->dimension();
   
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];
      
	    Transformation(T,iq);
	    GetFem()->ReInit(T);
	    if(dim==2)
	      {
		Vertex1d xi;
		xi.x() = 0.5;
		GetFem()->point_boundary(ile,xi);
	      }
	    else
	      {
		Vertex2d xi;
		xi.x() = 0.5;
		xi.y() = 0.5;
		GetFem()->point_boundary(ile,xi);
	      }
	    if(a[iq] == 1.)
	      a[iq] = GetFem()->G();
	    else
	      a[iq] += GetFem()->G();
	  }
      }
  }
  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::GetMassDiag(DoubleVector& a) const
  {
    a.resize(GetMesh()->nnodes());
    nmatrix<double> T;
    DoubleVector F;
    for(int iq=0;iq<GetMesh()->ncells();++iq)
      {
	Transformation(T,iq);
	GetFem()->ReInit(T);

	GetIntegrator()->IntegrateMassDiag(F,*GetFem());
	
	//Auf Globalen Vektor verteielen
	IntVector indices = GetLocalIndices(iq);
	for(int ii=0; ii<indices.size(); ii++) 
	  {
	    int i = indices[ii];
	    a[i] += F[ii];
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::GetBoundaryMassDiag(DoubleVector& a) const
  {
    a.resize(GetMesh()->nnodes());
    a.equ(1.);
    
    const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
    assert(GMP);
    std::set<int> Colors =  GMP->GetColors();
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& bv = *GMP->VertexOnBoundary(col);
	for(int i=0;i<bv.size();i++)
	  {
	    a[bv[i]] = 0;
	  }  
      }
    
    DoubleVector F;
    nmatrix<double> T;
    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	const IntVector& q = *GetMesh()->CellOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalOnBoundary(col);
	for (int i=0; i<q.size(); i++)
	  {
	    int iq  = q[i];
	    int ile = l[i];

	    Transformation(T,iq);
	    GetFem()->ReInit(T);

	    GetIntegrator()->IntegrateBoundaryMassDiag(F,*GetFem(),ile,col);

	    //Auf Globalen Vektor verteielen
	    IntVector indices = GetLocalIndices(iq);
	    for(int ii=0; ii<indices.size(); ii++) 
	      {
		int i = indices[ii];
		a[i] += F[ii];
	      }
	  }
      }
  }

  template<int DIM>
  void Q1<DIM>::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
  {
    IntVector indices = GetLocalIndices(iq);
    HN->CondenseHanging(E,indices);
    IntVector::const_iterator  start = indices.begin();
    IntVector::const_iterator  stop  = indices.end();
    A.entry(start,stop,__E,s);
  }
  
  template<int DIM>
  void Q1<DIM>::StrongDirichletMatrix(MatrixInterface& A, int col, const vector<int>& comp) const
  {
    const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
    assert(GMP);
    const IntVector& bv = *GMP->VertexOnBoundary(col);
    for(int i=0;i<bv.size();i++)
      {
	A.dirichlet(bv[i], comp);
      }  
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const vector<int>& comp) const
  {
    const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
    assert(GMP);
    const IntVector& bv = *GMP->VertexOnBoundary(col);
    for(int i=0;i<bv.size();i++)
      {
	A.dirichlet_only_row(bv[i], comp);
      }  
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::StrongDirichletVectorZero(GlobalVector& u, int col, const vector<int>& comp) const
  {
    const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
    assert(GMP);
    const IntVector& bv = *GMP->VertexOnBoundary(col);

    for(int ii=0;ii<comp.size();ii++)
      {
	int c = comp[ii];
	if(c<0) {
	  cerr << "negative component: " << c << endl;
	  abort();
	} else if(c>=u.ncomp()){
	  cerr << "unknown component: " << c << endl;
	  abort();
	}
      }

    for(int i=0;i<bv.size();i++)
      {
	int index = bv[i];
	for(int ii=0;ii<comp.size();ii++)
	  {
	    u( index,comp[ii] ) = 0.;
	  }
      }  
  }


  template<>
  void Q1<2>::InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
  {
    const IntVector& vo2n = *GetMesh()->Vertexo2n();
    nvector<bool> habschon(GetMesh()->nnodes(),0);  
    
    assert(vo2n.size()==uold.n());
    assert(GetMesh()->nnodes()==u.n());
    assert(u.ncomp()==uold.ncomp());
    
    for(int i=0;i<vo2n.size();i++)
      {
	int in = vo2n[i];
	
	if(in>=0) 
	  {
	    u.equ_node(in,1.,i,uold);
	    habschon[in] = 1;
	  }
      }
    nvector<fixarray<3,int> > nodes(4);
    nodes[0][0] = 1; nodes[0][1] = 0;  nodes[0][2] = 2;
    nodes[1][0] = 3; nodes[1][1] = 0;  nodes[1][2] = 6;
    nodes[2][0] = 5; nodes[2][1] = 2;  nodes[2][2] = 8;
    nodes[3][0] = 7; nodes[3][1] = 6;  nodes[3][2] = 8;
    
    const PatchMesh* PM = dynamic_cast<const PatchMesh*>(GetMesh());
    assert(PM);
    
    for(int iq=0;iq<PM->npatches();++iq)
      {
	IntVector vi =* PM->IndicesOfPatch(iq);
	
	for(int j=0; j<nodes.size(); j++)
	  {
	    int v  = vi[nodes[j][0]];
	    int v1 = vi[nodes[j][1]];
	    int v2 = vi[nodes[j][2]];
	    assert(habschon[v1]);
	    assert(habschon[v2]);
          if (habschon[v]==0) 
            {
              u.equ_node(v,0.5,v1,uold);
              u.add_node(v,0.5,v2,uold);
              habschon[v] = 1;
            }
	  }
	int v = vi[4];
	if (habschon[v]==0)
	  {
	    u.equ_node(v,0.25,vi[0],uold);
	    u.add_node(v,0.25,vi[2],uold);	  
	    u.add_node(v,0.25,vi[6],uold);	  
	    u.add_node(v,0.25,vi[8],uold);	  
	    habschon[v] = 1;
	  }
      }
  }
  
  template<>
  void Q1<3>::InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
  {
    const IntVector& vo2n = *GetMesh()->Vertexo2n();
    nvector<bool> habschon(GetMesh()->nnodes(),0);  
    
    assert(vo2n.size()==uold.n());
    
    for(int i=0;i<vo2n.size();i++)
      {
	int in = vo2n[i];
	
	if(in>=0) 
	  {
	    u.equ_node(in,1.,i,uold);
	    habschon[in] = 1;
	  }
      }
    nvector<fixarray<3,int> > nodes(12);
    nodes[0][0] = 1;    nodes[0][1] = 0;     nodes[0][2] = 2;
    nodes[1][0] = 3;    nodes[1][1] = 0;     nodes[1][2] = 6;
    nodes[2][0] = 5;    nodes[2][1] = 2;     nodes[2][2] = 8;
    nodes[3][0] = 7;    nodes[3][1] = 6;     nodes[3][2] = 8;
    nodes[4][0] = 1+18; nodes[4][1] = 0+18;  nodes[4][2] = 2+18;
    nodes[5][0] = 3+18; nodes[5][1] = 0+18;  nodes[5][2] = 6+18;
    nodes[6][0] = 5+18; nodes[6][1] = 2+18;  nodes[6][2] = 8+18;
    nodes[7][0] = 7+18; nodes[7][1] = 6+18;  nodes[7][2] = 8+18;
    nodes[8][0] = 9;    nodes[8][1] = 0;     nodes[8][2] = 18;
    nodes[9][0] = 11;   nodes[9][1] = 2;     nodes[9][2] = 20;
    nodes[10][0] = 15;  nodes[10][1] = 6;    nodes[10][2] = 24;
    nodes[11][0] = 17;  nodes[11][1] = 8;    nodes[11][2] = 26;
    
    
    nvector<fixarray<5,int> > w(6);
    w[0][0] = 4;  w[0][1] = 0;  w[0][2] = 2;  w[0][3] = 6;  w[0][4] = 8;
    w[1][0] = 12; w[1][1] = 0;  w[1][2] = 18; w[1][3] = 6;  w[1][4] = 24;
    w[2][0] = 14; w[2][1] = 2;  w[2][2] = 8;  w[2][3] = 20; w[2][4] = 26;
    w[3][0] = 16; w[3][1] = 6;  w[3][2] = 8;  w[3][3] = 24; w[3][4] = 26;
    w[4][0] = 10; w[4][1] = 0;  w[4][2] = 2;  w[4][3] = 18; w[4][4] = 20;
    w[5][0] = 22; w[5][1] = 18; w[5][2] = 20; w[5][3] = 24; w[5][4] = 26;
    
    const PatchMesh* PM = dynamic_cast<const PatchMesh*>(GetMesh());
    assert(PM);
    
    for(int iq=0;iq<PM->npatches();++iq)
      {
	IntVector vi = *PM->IndicesOfPatch(iq);
	
	for(int j=0; j<nodes.size(); j++)
	  {
	    int v  = vi[nodes[j][0]];
	    int v1 = vi[nodes[j][1]];
	    int v2 = vi[nodes[j][2]];
	    assert(habschon[v1]);
	    assert(habschon[v2]);
	    if (habschon[v]==0) 
	      {
		u.equ_node(v,0.5,v1,uold);
		u.add_node(v,0.5,v2,uold);
		habschon[v] = 1;
	      }
	  }
	for(int j=0; j<w.size(); j++)
	  {
	    int v  = vi[w[j][0]];
	    int v1 = vi[w[j][1]];
	    int v2 = vi[w[j][2]];
	    int v3 = vi[w[j][3]];
	    int v4 = vi[w[j][4]];
	    assert(habschon[v1]);
	    assert(habschon[v2]);
	    assert(habschon[v3]);
	    assert(habschon[v4]);
	    if (habschon[v]==0) 
	      {
		u.equ_node(v,0.25,v1,uold);
		u.add_node(v,0.25,v2,uold);
		u.add_node(v,0.25,v3,uold);
		u.add_node(v,0.25,v4,uold);
		habschon[v] = 1;
	      }
	  }
	int v = vi[13];
	if (habschon[v]==0)
	  {
	    u.equ_node(v,0.125,vi[0],uold);
	    u.add_node(v,0.125,vi[2],uold);	  
	    u.add_node(v,0.125,vi[6],uold);	  
	    u.add_node(v,0.125,vi[8],uold);	  
	    u.add_node(v,0.125,vi[18],uold);
	    u.add_node(v,0.125,vi[20],uold);	  
	    u.add_node(v,0.125,vi[24],uold);	  
	    u.add_node(v,0.125,vi[26],uold);	  
	    habschon[v] = 1;
	  }
      }
  }
  
  /*-----------------------------------------*/
  
  template<int DIM>
  void Q1<DIM>::HNAverage(GlobalVector& x) const
  {
    HN->Average(x);
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::HNDistribute(GlobalVector& x) const
  {
    HN->Distribute(x);
  }

  /* ----------------------------------------- */

  template<int DIM>
  void Q1<DIM>::HNZero(GlobalVector& x) const
  {
    HN->Zero(x);
  }

  /* ----------------------------------------- */

  template<int DIM>
  bool Q1<DIM>::HNZeroCheck(const GlobalVector& x) const
  {
    return HN->ZeroCheck(x);
  }


  /* ----------------------------------------- */
  
  template<>
  HNStructureInterface* Q1<2>::NewHNStructure()
  {
    return new HNStructureQ12d;
  }
  /* ----------------------------------------- */
  
  template<>
  HNStructureInterface* Q1<3>::NewHNStructure()
  {
    return new HNStructureQ13d;
  }
  
  template<int DIM>
  void Q1<DIM>::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT) 
  {
    {
      MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
      if(IP)
	{
	  IP->BasicInit(MT);
	  return;
	}
    }
    
    MgInterpolatorMatrix* IP = dynamic_cast<MgInterpolatorMatrix*>(I);
    assert(IP);
    const GascoigneMeshTransfer* GT = dynamic_cast<const GascoigneMeshTransfer*>(MT);
    assert(GT);
    
    const map<int,fixarray<2,int> >& zweier = GT->GetZweier();
    const map<int,fixarray<4,int> >& vierer = GT->GetVierer();
    const map<int,fixarray<8,int> >& achter = GT->GetAchter();
    const IntVector& c2f    = GT->GetC2f();
    
    int n  = c2f.size() +   zweier.size() +   vierer.size() +   achter.size();
    int nt = c2f.size() + 2*zweier.size() + 4*vierer.size() + 8*achter.size();
    
    ColumnStencil& ST = IP->GetStencil();
    DoubleVector& val = IP->GetAlpha();
    
    SparseStructure SS;
    
    SS.build_begin(n);
    for(int i=0;i<c2f.size();i++)
      {
	assert(c2f[i]>=0);
	
	SS.build_add(c2f[i],i);
      }
    for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
	p!=zweier.end();p++) 
      {
	int il = p->first;
	fixarray<2,int> n2 = p->second;
	for(int ii=0;ii<2;ii++) SS.build_add(il,n2[ii]);
      }
    for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
	p!=vierer.end();p++) {
      int il = p->first;
      fixarray<4,int> n4 = p->second;
      for(int ii=0;ii<4;ii++) SS.build_add(il,n4[ii]);
    }
    for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
	p!=achter.end();p++) {
      int il = p->first;
      fixarray<8,int> n8 = p->second;
      for(int ii=0;ii<8;ii++) SS.build_add(il,n8[ii]);
    }
    SS.build_end();
    
    assert(nt==SS.ntotal());
    
    ST.memory(&SS);
    
    val.reservesize(nt);
    
    for(int i=0;i<c2f.size();i++)
      {
	// ich weiss nicht, ob das richtig ist !!!!!
	int pos = ST.Find(c2f[i],i);
	assert(pos>=0);
	val[pos] = 1.;
      }
    for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
	p!=zweier.end();p++) 
      {
	int il = p->first;
	fixarray<2,int> n2 = p->second;
	val[ST.Find(il,n2[0])] = 0.5;
	val[ST.Find(il,n2[1])] = 0.5;
      }
    for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
	p!=vierer.end();p++) {
      int il = p->first;
      fixarray<4,int> n4 = p->second;
      val[ST.Find(il,n4[0])] = 0.25;
      val[ST.Find(il,n4[1])] = 0.25;
      val[ST.Find(il,n4[2])] = 0.25;
      val[ST.Find(il,n4[3])] = 0.25;
    }
    for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
	p!=achter.end();p++) {
      int il = p->first;
      fixarray<8,int> n8 = p->second;
      for (int i=0; i<8; i++)
	{
	  val[ST.Find(il,n8[i])] = 0.125;
	}
    }
  }


  /* ----------------------------------------- */

  template<>
  void Q1<2>::StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const vector<int>& comp, double d) const 
  {
    const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
    assert(GMP);
    DoubleVector ff(u.ncomp(),0.);
    const IntVector& bv = *GMP->VertexOnBoundary(col);
    
    FemData QH;
    
    GlobalToGlobalData();
    BF.SetParameterData(__QP);
    
    for(int ii=0;ii<comp.size();ii++)
      {
	int c = comp[ii];
	if(c<0) {
	  cerr << "negative component: " << c << endl;
	  abort();
	} else if(c>=u.ncomp()){
	  cerr << "unknown component: " << c << endl;
	  abort();
	}
      }
    
    for(int i=0;i<bv.size();i++)
      {
	int index = bv[i];
	
	QH.clear();
	GlobalData::const_iterator p=GetDataContainer().GetNodeData().begin();
	for(; p!=GetDataContainer().GetNodeData().end(); p++)
	  {
	    QH[p->first].resize(p->second->ncomp());
	    for(int c=0; c<p->second->ncomp(); c++)
	      {
		QH[p->first][c].m() = p->second->operator()(index,c);
	      }
	  }
	
	BF.SetFemData(QH);
	
	const Vertex2d& v = GMP->vertex2d(index);
	
	BF(ff,v,col);
	for(int iii=0;iii<comp.size();iii++)
	  {
	    int c = comp[iii];
	    u(index,c) = d * ff[c];
	  }
      }
  }

  template<>
  void Q1<3>::StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const vector<int>& comp, double d) const 
  {
    const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
    assert(GMP);
    DoubleVector ff(u.ncomp(),0.);
    const IntVector& bv = *GMP->VertexOnBoundary(col);
    
    FemData QH;
    
    GlobalToGlobalData();
    BF.SetParameterData(__QP);
    
    for(int ii=0;ii<comp.size();ii++)
      {
	int c = comp[ii];
	if(c<0) {
	  cerr << "negative component: " << c << endl;
	  abort();
	} else if(c>=u.ncomp()){
	  cerr << "unknown component: " << c << endl;
	  abort();
	}
      }
    
    for(int i=0;i<bv.size();i++)
      {
	int index = bv[i];
	
	QH.clear();
	GlobalData::const_iterator p=GetDataContainer().GetNodeData().begin();
	for(; p!=GetDataContainer().GetNodeData().end(); p++)
	  {
	    QH[p->first].resize(p->second->ncomp());
	    for(int c=0; c<p->second->ncomp(); c++)
	      {
		QH[p->first][c].m() = p->second->operator()(index,c);
	      }
	  }
	
	BF.SetFemData(QH);
	
	const Vertex3d& v = GMP->vertex3d(index);
	
	BF(ff,v,col);
	for(int iii=0;iii<comp.size();iii++)
	  {
	    int c = comp[iii];
	    u(index,c) = d * ff[c];
	  }
      }
  }
  
  
  
  

  template class Q1<2>;
  template class Q1<3>;
  
}
