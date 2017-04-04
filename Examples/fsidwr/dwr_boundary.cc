
#include  "dwr_boundary.h"
#include  "pi.h"
#include  "dwrfem.h"
#include  "alediscretization.h"
#include <iostream>
#include "stdsolver.h"

using namespace std;
namespace Gascoigne
{

  void DWRBoundary::reinit_element_2d(int en, const nvector<int>& indices, 
				    HASHMAP<int, std::vector<int> >& solid_interface_cells, 
				    HASHMAP<int, std::vector<int> >& fluid_interface_cells,
				    HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
				    HASHSET<int> & interface_nodes,
				    set<int>& fluid_nodes, set<int>& solid_nodes)
  {
    Chi chi;

    std::string __solid_type;
    DataFormatHandler DFH;
    DFH.insert("solid_type",&__solid_type);
    FileScanner FS(DFH,_paramfile,"Equation");

    chi.BasicInit(__solid_type);
    
    
    const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*> (__S->GetMesh());
    assert(M);

    int nf=0;
    int ns=0;
    
    vector<int> ni;
    for (int i=0;i<indices.size();++i)
      {
	int domain = chi(M->vertex2d(indices[i]));
	if (domain>0)  
	  {
	    ++ns;
	    solid_nodes.insert(indices[i]);
	  }
	if (domain<0) 
	  {
	    ++nf;
	    fluid_nodes.insert(indices[i]);
	  }
	if (domain==0) 
	  {
	    ni.push_back(i);
	    fluid_nodes.insert(indices[i]);
	    solid_nodes.insert(indices[i]);
	    interface_nodes.insert(indices[i]);
	  }
      }
    
    if ((ns>0)&&(nf>0))
      {

	cerr << "Geht nicht, fluid & solid!" << endl;

	for (int i=0;i<indices.size();++i)
	  cerr << M->vertex2d(indices[i])<< "\t" << chi(M->vertex2d(indices[i])) << endl;

	abort();
      }
    if (ni.size()>0)
      {
	if (ns>0)      
	  {
	    solid_interface_cells[en]=ni;
	    solid_cells.insert(en);
	  }
	else if (nf>0) 
	  {
	    fluid_interface_cells[en]=ni;
	    fluid_cells.insert(en);
	  }
	else 
	  {
	    solid_interface_cells[en]=ni;
	    cout << "Element has interface everywhere!!!" << endl;
	  }

	//	if (nf>0) cout << indices << "\t\t" << ni << endl;
	
      }
    else
      {
	if (ns==indices.size()) solid_cells.insert(en);
	else if (nf==indices.size()) fluid_cells.insert(en);
	else abort();
      }
  }

  void DWRBoundary::reinit_element_3d(int en, const nvector<int>& indices, 
				    HASHMAP<int, std::vector<int> >& solid_interface_cells, 
				    HASHMAP<int, std::vector<int> >& fluid_interface_cells,
				    HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
				    HASHSET<int> & interface_nodes,
				    set<int>& fluid_nodes, set<int>& solid_nodes)
  {
    Chi chi;

    std::string __solid_type;
    DataFormatHandler DFH;
    DFH.insert("solid_type",&__solid_type);
    FileScanner FS(DFH,_paramfile,"Equation");

    chi.BasicInit(__solid_type);
    
    
    const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (__S->GetMesh());
    assert(M);

    int nf=0;
    int ns=0;
    
    vector<int> ni;
    for (int i=0;i<indices.size();++i)
      {
	int domain = chi(M->vertex3d(indices[i]));
	if (domain>0)  
	  {
	    ++ns;
	    solid_nodes.insert(indices[i]);
	  }
	if (domain<0) 
	  {
	    ++nf;
	    fluid_nodes.insert(indices[i]);
	  }
	if (domain==0) 
	  {
	    ni.push_back(i);
	    fluid_nodes.insert(indices[i]);
	    solid_nodes.insert(indices[i]);
	    interface_nodes.insert(indices[i]);
	  }
      }
    
    if ((ns>0)&&(nf>0))
      {
	cerr << "Geht nicht, fluid & solid!" << endl;
	abort();
      }
    if (ni.size()>0)
      {
	if (ns>0)      
	  {
	    solid_interface_cells[en]=ni;
	    solid_cells.insert(en);
	  }
	else if (nf>0) 
	  {
	    fluid_interface_cells[en]=ni;
	    fluid_cells.insert(en);
	  }
	else 
	  {
	    solid_interface_cells[en]=ni;
	    cout << "Element has interface everywhere!!!" << endl;
	  }

	//	if (nf>0) cout << indices << "\t\t" << ni << endl;
	
      }
    else
      {
	if      (ns==indices.size()) solid_cells.insert(en);
	else if (nf==indices.size()) fluid_cells.insert(en);
	else 
	  {
	    cerr << ns << " " << nf << "\t" << indices.size() << endl;
	    
	    abort();
	  }
	
      }
  }
  

  void DWRBoundary::ReInitInterface(DiscretizationInterface* DISC)
  {
    AleBaseDiscretization* ALEDISC = dynamic_cast<AleBaseDiscretization*> (DISC);
    assert(ALEDISC);
    HASHMAP<int, std::vector<int> >& solid_interface_cells = ALEDISC->GetSolidInterfaceCells();
    HASHMAP<int, std::vector<int> >& fluid_interface_cells = ALEDISC->GetFluidInterfaceCells();
    HASHSET<int>&                    interface_nodes       = ALEDISC->GetInterfaceNodes();
    HASHSET<int>&                    fluid_cells           = ALEDISC->GetFluidCells();
    HASHSET<int>&                    solid_cells           = ALEDISC->GetSolidCells();
    vector<int>&                     fluid_l2g             = ALEDISC->GetFluidL2G();
    vector<int>&                     solid_l2g             = ALEDISC->GetSolidL2G();
    HASHMAP<int,int>&                fluid_g2l             = ALEDISC->GetFluidG2L();
    HASHMAP<int,int>&                solid_g2l             = ALEDISC->GetSolidG2L();

    set<int> fluid_nodes, solid_nodes;


    solid_interface_cells.clear();
    fluid_interface_cells.clear();
    interface_nodes.clear();
    fluid_cells.clear();
    solid_cells.clear();
    
    
    int dim = __S->GetMesh()->dimension();

    if (dim==2)
      {	
	const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*> (__S->GetMesh());
	assert(M);
	
	if ((DISC->GetName()=="Q1 Ale 2d Lps")||
	    (DISC->GetName()=="Q1 Ale 2d"))
	  for (int c=0;c<M->ncells();++c)
	    reinit_element_2d(c, M->IndicesOfCell(c), solid_interface_cells, 
			      fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	
	else if ((DISC->GetName()=="Q2 Ale Lps 2d")||
		 (DISC->GetName()=="Q2 Ale 2d"))
	  for (int c=0;c<M->npatches();++c)
	    reinit_element_2d(c, *(M->IndicesOfPatch(c)), solid_interface_cells, 
			      fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	else abort();

      }
    else if (dim==3)
      {	
	const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (__S->GetMesh());
	assert(M);
	
	if (DISC->GetName()=="Q1 Ale Lps 3d")
	  for (int c=0;c<M->ncells();++c)
	    reinit_element_3d(c, M->IndicesOfCell(c), solid_interface_cells, fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	else if (DISC->GetName()=="Q2 Ale Lps 3d")
	  for (int c=0;c<M->npatches();++c)
	    reinit_element_3d(c, *(M->IndicesOfPatch(c)), solid_interface_cells, fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	else abort();

      }
    else abort();

    // Nodes Fluid & Solid,  local <-> global (fluid nodes include interface and same for solid)
    fluid_l2g.clear();
    solid_l2g.clear();
    fluid_g2l.clear();
    solid_g2l.clear();

    // l2g
    for (set<int>::const_iterator it = fluid_nodes.begin();it!=fluid_nodes.end();++it)
      fluid_l2g.push_back(*it);
    for (set<int>::const_iterator it = solid_nodes.begin();it!=solid_nodes.end();++it)
      solid_l2g.push_back(*it);

    // g2l
    for (int i=0;i<fluid_l2g.size();++i) fluid_g2l[fluid_l2g[i]] = i;
    for (int i=0;i<solid_l2g.size();++i) solid_g2l[solid_l2g[i]] = i;


  }
























  /*--------------------------------------------------------*/
  
  DWRBoundary::DWRBoundary(MultiLevelSolverInterface* MS, 
			   const std::string primalproblem,  const std::string dualproblem,const ParamFile* pf) 
    :  _paramfile(pf),  __MS(dynamic_cast<StdMultiLevelSolver*> (MS)), __primalproblem(primalproblem), __dualproblem(dualproblem)
  { 
    assert(__MS);
    __S = __MS->GetSolver();
    assert(__S);
    __D = __S->GetDiscretization();
    assert(__D);

    assert(__MS->GetProblemContainer()->GetProblem(__primalproblem));
    assert(__MS->GetProblemContainer()->GetProblem(__dualproblem));

    //__higher = new DWR_Discretization("Q1patch","Q1patch", "Q2","Q2");
    //__higher = new DWR_Discretization("Q2",     "Q1patch", "Q2",     "Q2");

    
    __lower  = new DWR_Discretization("Q1patch","Q1patch", "Q1patch","Q1patch");

    
    //    __higher = new AleQ22d;    
    abort();
    
    __higher->BasicInit(__S->GetParamfile());
    __higher->ReInit(__S->GetMesh());
    __higher->SetDataContainer(__S->GetDiscretization()->GetDataContainer());
    ReInitInterface(__higher);

    
    __lower->BasicInit(__S->GetParamfile());
    __lower ->ReInit(__S->GetMesh());
    __lower ->SetDataContainer(__S->GetDiscretization()->GetDataContainer());
  }

  
  /*--------------------------------------------------------*/
  
  double DWRBoundary::ScalarProduct(DoubleVector& eta, const GlobalVector& f,  const GlobalVector& z) const
  {
    double res = 0.0;
    for(int i=0; i<z.n(); i++)
      for (int c=0; c<z.ncomp(); c++)
	{
	  res += f(i,c) * z(i,c);
	  eta[i] += (f(i,c)*z(i,c));
	}
       
    return res;
  }
  
  /*--------------------------------------------------------*/
  
  double DWRBoundary::ScalarProduct(DoubleVector& eta, const VectorInterface& gf, const VectorInterface& gz) const
  {
    const GlobalVector& f = __S->GetGV(gf);
    const GlobalVector& z = __S->GetGV(gz);
    
    return ScalarProduct(eta,f,z);
  }
  
  /*--------------------------------------------------------*/
  
  double DWRBoundary::ScalarProductWithFluctuations(DoubleVector& eta, const VectorInterface& gf, const VectorInterface& gz) const
  {
    double X =  ScalarProduct(eta,gf,gz);

    const GlobalVector& f = __S->GetGV(gf);
    const GlobalVector& z = __S->GetGV(gz);
    
    GlobalVector dz(f.ncomp(),f.n());
    
    dz.zero();
    Pi pi;
    pi.Init(__S->GetMesh());
    pi.vmult(dz,z);

    eta.zero();
    ScalarProduct(eta,f,dz);
    //    ScalarProduct(eta,f,z);
    
    return X;
  }

  /*-------------------------------------------------------*/
  
  void DWRBoundary::ResidualsHigher(VectorInterface& f, const VectorInterface& u,DiscretizationInterface* D, double s)
  {
    assert(__D);
    assert(D);
    
    // if (D!=__D)
    //   dynamic_cast<const DWR_Discretization* >(D)->HNAverageAnsatz(const_cast<GlobalVector& >(__S->GetGV(u)));
    // else
    //   D->HNAverage(const_cast<GlobalVector&>(__S->GetGV(u)));
	
    if (D!=__D)
      __S->SetDiscretization(*D);
    

    __S->Rhs (f,   -s);
    __S->Form(f,u,  s);



    // if (D!=__D)
    //   dynamic_cast<const DWR_Discretization*>(D)->HNDistributeTest(const_cast<GlobalVector&>(__S->GetGV(f)));
    // else
    //   D->HNDistribute(const_cast<GlobalVector&>(__S->GetGV(f)));
    
    if (D!=__D)
      __S->SetDiscretization(*__D);
    
    // if (D!=__D)
    //   dynamic_cast<const DWR_Discretization*>(D)->HNZeroAnsatz(const_cast<GlobalVector&>(__S->GetGV(u)));
    // else 
    //   D->HNZero(const_cast<GlobalVector&>(__S->GetGV(u)));
  }
  
  /*--------------------------------------------------------*/
  
  double DWRBoundary::Estimator(DoubleVector& eta, VectorInterface& f, VectorInterface& u, VectorInterface& z)
  {
    eta.resize(__S->GetGV(z).n());
    eta.zero();
    

    __S->SetProblem(*__MS->GetProblemContainer()->GetProblem(__dualproblem));
    __higher->AddNodeVector("U",&__S->GetGV(u));
    __D->AddNodeVector("U",&__S->GetGV(u));

    __S->Zero(f);
    ResidualsHigher(f,z,__higher,-0.5);
    ResidualsHigher(f,z,__D,  0.5);

    __higher->DeleteNodeVector("U");
    __D->DeleteNodeVector("U");
    double rhostern =  ScalarProductWithFluctuations(eta,f,u);

    __S->SetProblem(*__MS->GetProblemContainer()->GetProblem(__primalproblem));

    __S->Zero(f);
    ResidualsHigher(f,u,__higher,-0.5);
    ResidualsHigher(f,u,__D,  0.5);

    double rho      =  ScalarProductWithFluctuations(eta,f,z);

    cout << "Estimator: primal residuals: " << rho << "\t dual residuals: " << rhostern << endl;
    cout << "Full Estimator: " << 0.5 * rho + 0.5 * rhostern << endl;
    

    return rho + rhostern;
  }
  
  /*--------------------------------------------------------*/
  
  double DWRBoundary::EstimatorEasy(DoubleVector& eta, VectorInterface& f, VectorInterface& u, VectorInterface& z) 
  {
    eta.resize(__S->GetGV(z).n());
    eta.zero();
    
    
    __S->Zero(f);
    ResidualsHigher(f,u,__higher,-0.5);
    ResidualsHigher(f,u,__D,  0.5);
    
    
    double rho      =  ScalarProductWithFluctuations(eta,f,z);
    
    return 2.0 * rho;
  }

  
  /*--------------------------------------------------------*/
  
  double DWRBoundary::EstimatorAdjoint(DoubleVector& eta, VectorInterface& f, VectorInterface& u, VectorInterface& z) 
  {
    eta.resize(__S->GetGV(z).n());
    eta.zero();
    
    __S->SetProblem(*__MS->GetProblemContainer()->GetProblem(__dualproblem));
    __S->Zero(f);

    __higher->AddNodeVector("U",&__S->GetGV(u));
    __D->AddNodeVector("U",&__S->GetGV(u));
    __ADJOINT = true;
    ResidualsHigher(f,z,__higher,-0.5);
    ResidualsHigher(f,z,__D,      0.5);
    __ADJOINT = false;
    __higher->DeleteNodeVector("U");
    __D->DeleteNodeVector("U");

    //    __S->SetBoundaryVectorZero(f);

    __S->SetProblem(*__MS->GetProblemContainer()->GetProblem(__primalproblem));
    double rho      =  ScalarProductWithFluctuations(eta,f,u);
    
    return 2.0 * rho;
  }
  
  
  
/*--------------------------------------------------------*/
}
