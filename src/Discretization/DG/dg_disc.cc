#include "dg_disc.h"
#include  "sparsestructure.h"
#include <fstream>

using namespace std;

namespace Gascoigne
{

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////  DISC
  //////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////// DOF-Handling

  template<int DIM>
   void Disc<DIM>::Structure(SparseStructureInterface* SI) const
  {
    assert(_dof_handler);
    SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
    assert(S);
    
    S->build_begin(GetDofHandler()->ndofs());

    // first, inner-element-coupling
    vector<int> indices(GetDofHandler()->ndofs_per_element());
    for (int el=0;el<GetDofHandler()->nelements();++el)
      {
	GetDofHandler()->GetIndices(el,indices);
	S->build_add(indices.begin(), indices.end());
      }
    // edges... // hd?
    
    S->build_end();  

  }

  
  template<int DIM>
  void Disc<DIM>::Transformation(FemInterface::Matrix& T, int iq) const
  {
    IntVector indices = GetMesh()->IndicesOfCell(iq);
    T.memory(DIM,indices.size());

    for(int ii=0;ii<indices.size();ii++)
      {
	const Vertex<DIM>& v = GetMesh()->vertex(indices[ii]);
	for (int d=0;d<DIM;++d)
	  T(d,ii) = v[d];
      }
  }
 
  ////////////////////////////////////////////////// Integration
  
  template<int DIM>
  void Disc<DIM>::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
  {
    assert(_dof_handler);
    nmatrix<double> T;

#warning DGDisc::Form GlobalToGlobalData missing
    //    GlobalToGlobalData();
    //    EQ.SetParameterData(__QP);
    LocalData __QN,__QC;

    LocalVector U,F;
    GetDofHandler()->InitLocalVector(U, EQ.GetNcomp());
    GetDofHandler()->InitLocalVector(F, EQ.GetNcomp());

    FiniteElementType FEM;

    for(int el=0;el<GetDofHandler()->nelements();++el)
      {
	Transformation(T,el);
	FEM.ReInit(T);
	
	GetDofHandler()->GlobalToLocal(U,u,el);	  
	_integrator.Form(EQ,F,FEM,U,__QN,__QC);	  
	GetDofHandler()->LocalToGlobal(f,F,el,d); 
      }
  }
  
  template<int DIM>
  void Disc<DIM>::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
  {
    assert(_dof_handler);
    nmatrix<double> T;

#warning DGDisc::Form GlobalToGlobalData missing
    //    GlobalToGlobalData();
    //    RHS.SetParameterData(__QP);
    LocalData __QN,__QC;

    LocalVector F;
    GetDofHandler()->InitLocalVector(F, RHS.GetNcomp());

    FiniteElementType FEM;

    for(int el=0;el<GetDofHandler()->nelements();++el)
      {
	Transformation(T,el);
	FEM.ReInit(T);

#warning DGDisc::Rhs GlobalToGlobalData(eq) missing
	// GlobalToLocalData(eq);
	      
	_integrator.Rhs(RHS,F,FEM,__QN,__QC);	  
	GetDofHandler()->LocalToGlobal(f,F,el,s);
      }
  }


  
  template<int DIM>
  void Disc<DIM>::Matrix(MatrixInterface& A, const GlobalVector& u, const ProblemDescriptorInterface* PD, double d) const
  {
    assert(_dof_handler);

    nmatrix<double> T;
    const Equation& EQ = *(PD->GetEquation());

#warning DGDisc::Form GlobalToGlobalData missing
    //GlobalToGlobalData();
    //EQ.SetParameterData(__QP);
    LocalData __QN,__QC;
    
    LocalVector U;
    GetDofHandler()->InitLocalVector(U,EQ.GetNcomp());

    EntryMatrix E;
    
    FiniteElementType FEM;

    for(int el=0;el<GetDofHandler()->nelements();++el)
      {
	Transformation(T,el);
	FEM.ReInit(T);

	GetDofHandler()->GlobalToLocal(U,u,el);
	//EQ.cell(GetMesh(),iq,__U,__QN);
	_integrator.Matrix(EQ,E,FEM,U,__QN,__QC);
	GetDofHandler()->LocalToGlobal(A,E,el,d);
      }
#warning hanging nodes???
    //    HN->MatrixDiag(u.ncomp(),A);
  }



  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////  DG
  //////////////////////////////////////////////////////////////////////

  // add stuff for edges
  template<int DIM>
  void DGDisc<DIM>::Form(GlobalVector& f, const GlobalVector& u, const Equation& BEQ, double d) const
  {
    assert(dynamic_cast<const EdgeEquation*> (&BEQ));
    const EdgeEquation& EQ = dynamic_cast<const EdgeEquation&> (BEQ);
    Disc<DIM>::Form(f,u,BEQ,d);

    nmatrix<double> T1,T2;
    FiniteElementType FEM1,FEM2;


#warning DGDisc::Form GlobalToGlobalData missing
    //    GlobalToGlobalData();
    //    EQ.SetParameterData(__QP);
    LocalData __QN,__QC;
    
    LocalVector U1,U2,F1,F2;
    this->GetDofHandler()->InitLocalVector(U1, BEQ.GetNcomp());
    this->GetDofHandler()->InitLocalVector(U2, BEQ.GetNcomp());
    this->GetDofHandler()->InitLocalVector(F1, BEQ.GetNcomp());
    this->GetDofHandler()->InitLocalVector(F2, BEQ.GetNcomp());

    for (int ed=0;ed<this->GetDGDofHandler()->nedges();++ed)
      {
	const DGEdge& E = this->GetDGDofHandler()->edge(ed);
	int c1 = E.master(); assert(c1!=-1);
	int c2 = E.slave();

	this->Transformation(T1,c1);
	FEM1.ReInit(T1);
	this->GetDofHandler()->GlobalToLocal(U1,u,c1);

	if (c2!=-1)
	  {
	    this->Transformation(T2,c2);
	    FEM2.ReInit(T2);
	    this->GetDofHandler()->GlobalToLocal(U2,u,c2);
	  }

	_edge_integrator.EdgeForm(EQ,E, F1,F2,FEM1,FEM2,U1,U2,__QN,__QC);
	this->GetDofHandler()->LocalToGlobal(f,F1,c1,d);
	if (c2!=-1)
	  this->GetDofHandler()->LocalToGlobal(f,F2,c2,d); 
      }
  }



  
  
  // I/O
  template<int DIM>
  void DGDisc<DIM>::WriteVtk(std::string name, const GlobalVector& u) const
  {
    ofstream OUT(name.c_str());
    assert(OUT);
    assert(OUT.is_open());

    // write header
    OUT << "# vtk DataFile Version 2.0" << endl;
    OUT << "Gascoigne DG-Output," << endl;
    OUT << "ASCII" << endl;

    // Type of Elements
    OUT << "DATASET UNSTRUCTURED_GRID" << endl;

    //////// Write Mesh
    // nodes
    OUT << "POINTS " << this->ndofs() << " DOUBLE" << endl;
    for (int c=0;c<this->GetMesh()->ncells();++c)
      {
	vector<int> indices = this->GetMesh()->IndicesOfCell(c);
	for (int i=0;i<indices.size();++i)
	  {
	    OUT << this->GetMesh()->vertex(indices[i]);
	    if (DIM==2) OUT << " 0";
	    OUT << endl;
	  }
      }
    // elements
    OUT << "CELLS " << this->GetDofHandler()->nelements() << " " << (DIM==2?5:9)*this->GetDofHandler()->nelements() << endl;
    vector<int> ind(this->GetDofHandler()->ndofs_per_element());
    for (int c=0;c<this->GetDofHandler()->nelements();++c)
      {
	OUT << (DIM==2?4:8) << " ";
	this->GetDofHandler()->GetIndices(c,ind);
	OUT << ind << endl;
      }
    OUT << "CELL_TYPES " << this->GetMesh()->ncells() << endl;
    for (int c=0;c<this->GetMesh()->ncells();++c) 
      OUT << (DIM==2?8:11) << " ";
    OUT << endl;

    // Data
    assert(u.ncomp()==1);
    OUT << "POINT_DATA " << u.n()  << endl;
    OUT << "SCALARS u000 DOUBLE" << endl;
    OUT << "LOOKUP_TABLE default" << endl;
    for (int i=0;i<u.n();++i)
      OUT << u(i,0) << " ";
    
    
    // 
    OUT.close();
  }
  
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////  CG
  //////////////////////////////////////////////////////////////////////


  // add stuff for hanging nodes
  

  template class DGDisc<2>;
  template class DGDisc<3>;

  template class CGDisc<2>;
  template class CGDisc<3>;
}
