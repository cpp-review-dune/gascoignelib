#include "dwr_discretization.h" 
#include "dwr_integrator.h" 

using namespace std;

namespace Gascoigne
{
  DWR_Discretization::DWR_Discretization(string trans_ansatz, string base_ansatz,
					 string trans_trial,  string base_trial) : PatchDiscretization(), __FEM_ANSATZ(0), __FEM_TRIAL(0)
  {
    if      ( (trans_ansatz=="Q1")     &&(base_ansatz=="Q1") )      
      {
	__FEM_ANSATZ = new FiniteElement<2,1,TransQ1,      BaseQ12d>;
	HNAnsatz = new HNStructureQ12d;
      }
    else if ( (trans_ansatz=="Q2")     &&(base_ansatz=="Q2") )      
      {
	__FEM_ANSATZ = new FiniteElement<2,1,TransQ2,      BaseQ22d>;
	HNAnsatz = new HNStructureQ22d;
      }
    else if ( (trans_ansatz=="Q2")     &&(base_ansatz=="Q1patch") ) 
      {
	__FEM_ANSATZ = new FiniteElement<2,1,TransQ2,      BaseQ12dPatch>;
	HNAnsatz = new HNStructureQ12d;
      }
    else if ( (trans_ansatz=="Q1patch")&&(base_ansatz=="Q1patch") ) 
      {
	__FEM_ANSATZ = new FiniteElement<2,1,TransQ1Patch, BaseQ12dPatch>;
	HNAnsatz = new HNStructureQ12d;
      }
    
    else
      {
	cerr << "== not defined: " << trans_ansatz << "\t" << base_ansatz << endl;
	abort();
      }

    if      ( (trans_trial=="Q1")     &&(base_trial=="Q1") )      
      {
	__FEM_TRIAL = new FiniteElement<2,1,TransQ1,      BaseQ12d>;
	HNTest = new HNStructureQ12d;
      }
    else if ( (trans_trial=="Q1patch")&&(base_trial=="Q1patch") ) 
      {
	__FEM_TRIAL = new FiniteElement<2,1,TransQ1Patch, BaseQ12dPatch>;
	HNTest = new HNStructureQ12d;
      }
    else if ( (trans_trial=="Q1patch")&&(base_trial=="Q2") )      
      {
	__FEM_TRIAL = new FiniteElement<2,1,TransQ1Patch, BaseQ22d>;
	HNTest = new HNStructureQ22d;
      }
    else if ( (trans_trial=="Q2")     &&(base_trial=="Q2") )      
      {
	__FEM_TRIAL = new FiniteElement<2,1,TransQ2,      BaseQ22d>;
	HNTest = new HNStructureQ22d;
      }
	  
    else
      {
	cerr << "== not defined: " << trans_trial << "\t" << base_trial << endl;
	abort();
      }


  }
  
  void DWR_Discretization::HNAverageAnsatz   (GlobalVector& x) const
  {
    HNAnsatz->Average(x);
  }
  
  void DWR_Discretization::HNZeroAnsatz      (GlobalVector& x) const
  {
    HNAnsatz->Zero(x);
  }
  
  void DWR_Discretization::HNDistributeTest  (GlobalVector& x) const
  {
    HNTest->Distribute(x);
  }
  
  // --------------------------------------------------

  DWR_Discretization::~DWR_Discretization()
  {
    if (__FEM_ANSATZ) delete __FEM_ANSATZ;
    if (__FEM_TRIAL)  delete __FEM_TRIAL;

    delete HNAnsatz;
    delete HNTest;
  }

  // --------------------------------------------------
  
  void DWR_Discretization::ReInit(const MeshInterface* MP)
  {
    PatchDiscretization::ReInit(MP);
    // hanging nodes!!!!
    HNAnsatz->ReInit(MP);
    HNTest->ReInit(MP);
  }

  /* ----------------------------------------- */

  void DWR_Discretization::BasicInit(const ParamFile* pf)
  {
    assert(!GetIntegratorPointer());
    GetIntegratorPointer() = new DWR_Integrator<2>;
    GetIntegratorPointer()->BasicInit();

    // haengende knoten!

    // Finite ELemente erstellen???? Nein, im Constructor!
    
    // Das brauch ich gar nicht.
    //CellDiscretization::BasicInit(pf);
  }

  // --------------------------------------------------
  
  void DWR_Discretization::GlobalToLocalAnsatz(LocalVector& U, const GlobalVector& u, int ip) const
  {
    GlobalToLocalSingleAnsatz(U,u,ip);
    GlobalToLocalDataAnsatz(ip);
  }
  void DWR_Discretization::GlobalToLocalTrial(LocalVector& U, const GlobalVector& u, int ip) const
  {
    GlobalToLocalSingleTrial(U,u,ip);
    GlobalToLocalDataTrial(ip);
  }
  // --------------------------------------------------
  void DWR_Discretization::GlobalToLocalSingleAnsatz(LocalVector& U, const GlobalVector& u, int ip) const
  {
    nvector<int> indices = GetAnsatzIndices(ip);
    U.ReInit(u.ncomp(),indices.size());
    for(int ii=0; ii<indices.size(); ii++) 
      {
	int i = indices[ii];
	assert(i<u.n());
	assert(ii<U.n());
	U.equ_node(ii,i,u);
      }
  }
  void DWR_Discretization::GlobalToLocalSingleTrial(LocalVector& U, const GlobalVector& u, int ip) const
  {
    nvector<int> indices = GetTrialIndices(ip);
    U.ReInit(u.ncomp(),indices.size());
    for(int ii=0; ii<indices.size(); ii++) 
      {
	int i = indices[ii];
	assert(i<u.n());
	assert(ii<U.n());
	U.equ_node(ii,i,u);
      }
  }
  // --------------------------------------------------
  void DWR_Discretization::GlobalToLocalDataAnsatz(int ip) const
  {
    const GlobalData& gnd = GetDataContainer().GetNodeData();
    __QN.clear();
    GlobalData::const_iterator p=gnd.begin();
    for(; p!=gnd.end(); p++)
      GlobalToLocalSingleAnsatz(__QN[p->first],*p->second,ip);
    
    const GlobalData& gcd = GetDataContainer().GetCellData();
    assert(gcd.size()==0);
  }
  void DWR_Discretization::GlobalToLocalDataTrial(int ip) const
  {						
    const GlobalData& gnd = GetDataContainer().GetNodeData();
    __QN.clear();
    GlobalData::const_iterator p=gnd.begin();
    for(; p!=gnd.end(); p++)
      GlobalToLocalSingleTrial(__QN[p->first],*p->second,ip);
    
    const GlobalData& gcd = GetDataContainer().GetCellData();
    assert(gcd.size()==0);    
  }
  //
  void DWR_Discretization::LocalToGlobalAnsatz(GlobalVector& f, const LocalVector& F, int iq, double s) const
  {
    nvector<int> indices = GetAnsatzIndices(iq);
    for(int ii=0; ii<indices.size(); ii++) 
      {
	int i = indices[ii];
	assert(i<f.n());
	assert(ii<F.n());

	f.add_node(i,s,ii,F);
      }
  }
  void DWR_Discretization::LocalToGlobalTrial(GlobalVector& f, const LocalVector& F, int iq, double s) const
  {
    nvector<int> indices = GetTrialIndices(iq);
    for(int ii=0; ii<indices.size(); ii++) 
      {
	int i = indices[ii];
	assert(i<f.n());
	assert(ii<F.n());

	f.add_node(i,s,ii,F);
      }
  }
  //
  nvector<int> DWR_Discretization::GetAnsatzIndices(int iq) const 
  {
    assert(GetAnsatzFem());
    if (GetAnsatzFem()->n()==4)      return  GetMesh()->IndicesOfCell(iq);
    else if (GetAnsatzFem()->n()==9) return *GetPatchMesh()->IndicesOfPatch(iq);
    else abort();
  }
  nvector<int> DWR_Discretization::GetTrialIndices(int iq) const 
  {
    assert(GetTrialFem());
    if (GetTrialFem()->n()==4)      return  GetMesh()->IndicesOfCell(iq);
    else if (GetTrialFem()->n()==9) return *GetPatchMesh()->IndicesOfPatch(iq);
    else abort();
  }

  /* ----------------------------------------- */

  void DWR_Discretization::TransformationQ1(FemInterface::Matrix& T, int ip) const
  {
    int dim = GetMesh()->dimension(); assert(dim==2);
    int ne = GetMesh()->nodes_per_cell(ip);
    nvector<int> indices = GetPatchMesh()->CoarseIndices(ip);
    assert(ne==indices.size());
    assert(ne==4);
    T.memory(dim,ne);
    for(int ii=0;ii<ne;ii++)
      {
	Vertex2d v = GetMesh()->vertex2d(indices[ii]);
	T(0,ii) = v.x();               
	T(1,ii) = v.y();
      }
    
  }
  
  void DWR_Discretization::TransformationQ2(FemInterface::Matrix& T, int ip) const
  {
    int dim = GetPatchMesh()->dimension(); assert(dim==2);
    int ne = GetPatchMesh()->nodes_per_patch();
    assert(ne==9);
    nvector<int> indices = *GetPatchMesh()->IndicesOfPatch(ip);
    assert(ne==indices.size());
    
    T.memory(dim,ne);
    for(int ii=0;ii<ne;ii++)
      {
	Vertex2d v = GetPatchMesh()->vertex2d(indices[ii]);
	T(0,ii) = v.x();               
	T(1,ii) = v.y();
      }
  }

  void DWR_Discretization::InitElement(int ip, const FemInterface* F) const
  {
    nmatrix<double> T;
    if      (F->n_trafo()==4) TransformationQ1(T,ip);
    else if (F->n_trafo()==9) TransformationQ2(T,ip);
    else abort();

    F->ReInit(T);
  }

  // --------------------------------------------------

  void DWR_Discretization::Form(GlobalVector& f, const GlobalVector& u, 
				const Equation& EQ, double d) const
  {
    GlobalToGlobalData();
    EQ.SetParameterData(__QP);
    
    for(int ip=0;ip<GetPatchMesh()->npatches();++ip)
      {
	InitElement(ip, GetAnsatzFem());
	InitElement(ip, GetTrialFem());
	
	GlobalToLocalAnsatz(__U,u,ip);
       	GetDWRIntegrator()->Form(EQ,__F,*GetAnsatzFem(),*GetTrialFem(),__U,__QN,__QC);
	LocalToGlobalTrial(f,__F,ip,d);
      }
  }

  // --------------------------------------------------

  void DWR_Discretization::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const
  {
    GlobalToGlobalData();
    BE.SetParameterData(__QP);


    for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
      {
	int col = *p;
	
	const IntVector& q = *GetMesh()->PatchOnBoundary(col);
	const IntVector& l = *GetMesh()->LocalPatchOnBoundary(col);
	
	for(int i=0;i<q.size();++i)
	  {
	    int ip  = q[i];
	    int ile = l[i];

	    InitElement(ip, GetAnsatzFem());
	    InitElement(ip, GetTrialFem());
	
	    GlobalToLocalAnsatz(__U,u,ip);
	    GetDWRIntegrator()->BoundaryForm(BE,__F,*GetAnsatzFem(),*GetTrialFem(),__U,ile,col,__QN,__QC);
	    LocalToGlobalTrial(f,__F,ip,d);
	  }
      }
  }
  
  /*---------------------------------------------------*/
    

  void DWR_Discretization::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
  {
    GlobalToGlobalData();
    RHS.SetParameterData(__QP);

    for(int ip=0;ip<GetPatchMesh()->npatches();++ip)
      {
	InitElement(ip, GetAnsatzFem());
	InitElement(ip, GetTrialFem());
	
	GlobalToLocalDataAnsatz(ip);
	GetDWRIntegrator()->Rhs(RHS, __F, *GetAnsatzFem(), *GetTrialFem(), __QN, __QC);
	LocalToGlobalTrial(f,__F,ip,s);
      }
  }
  
}

#include "finiteelement.xx"

template class Gascoigne::FiniteElement<2,1,Gascoigne::Transformation2d<Gascoigne::BaseQ12dPatch>, Gascoigne::BaseQ12dPatch>;



