#include  "dwr.h"
#include  "pi-fsi.h"
#include  "alediscretization.h"
#include <set>
#include  "paramfile.h"
#include "stdsolver.h"
#include  "lpsequation.h"

extern bool __ESTIMATE;
extern bool __ADJOINT;

using namespace std;

namespace Gascoigne
{


  void Dwr::reinit_element_2d(int en, const nvector<int>& indices, 
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
    FileScanner FS(DFH,dynamic_cast<const StdSolver&> (S).GetParamfile(),"Equation");

    chi.BasicInit(__solid_type);
    
    
    const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*> (S.GetMesh());
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

  void Dwr::reinit_element_3d(int en, const nvector<int>& indices, 
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
    FileScanner FS(DFH,S.GetParamfile(),"Equation");

    chi.BasicInit(__solid_type);
    
    
    const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (S.GetMesh());
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
  

  void Dwr::ReInitInterface(DiscretizationInterface* DISC)
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
    
    
    int dim = S.GetMesh()->dimension();

    if (dim==2)
      {	
	const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*> (S.GetMesh());
	assert(M);
	
	if ((DISC->GetName()=="Q1 Ale 2d Lps")||
	    (DISC->GetName()=="Q1 Ale 2d"))
	  for (int c=0;c<M->ncells();++c)
	    reinit_element_2d(c, M->IndicesOfCell(c), solid_interface_cells, 
			      fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	
	else if ((DISC->GetName()=="Q2 Ale Lps 2d")||
		 (DISC->GetName()=="Q2 Ale 2d")||
		 (DISC->GetName()=="DwrFem Q1 Q2 Ale 2d"))
	  for (int c=0;c<M->npatches();++c)
	    reinit_element_2d(c, *(M->IndicesOfPatch(c)), solid_interface_cells, 
			      fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	else abort();

      }
    else if (dim==3)
      {	
	const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (S.GetMesh());
	assert(M);
	
	if (DISC->GetName()=="Q1 Ale Lps 3d")
	  for (int c=0;c<M->ncells();++c)
	    reinit_element_3d(c, M->IndicesOfCell(c), solid_interface_cells, fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	else if ((DISC->GetName()=="Q2 Ale Lps 3d")||
		 (DISC->GetName()=="DwrFem Q1 Q2 Ale 3d"))
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
  
  Dwr::Dwr(SolverInterface& SR, 
	   const ProblemDescriptorInterface* primal,
	   const ProblemDescriptorInterface* dual) :
    S(SR) , primalproblem(primal), dualproblem(dual)
  { 
    discretization = S.GetDiscretization();
  }
  
  /*--------------------------------------------------------*/

  double Dwr::ScalarProduct(DoubleVector& eta, const GlobalVector& f, 
			     const GlobalVector& z) const
  {
    for(int i=0; i<z.n(); i++)
      {
	for (int c=0; c<z.ncomp(); c++)
	  {
	    eta[i] += fabs(f(i,c)*z(i,c));
	  }
      } 
    return z * f;
  }
  
  /*--------------------------------------------------------*/
  
  double Dwr::ScalarProduct(DoubleVector& eta, const VectorInterface& gf, 
			    const VectorInterface& gz) const
  {
    const GlobalVector& f = S.GetGV(gf);
    const GlobalVector& z = S.GetGV(gz);
    
    return ScalarProduct(eta,f,z);
  }
  
  /*--------------------------------------------------------*/

  double Dwr::ScalarProductWithFluctuations(DoubleVector& eta, 
					    const VectorInterface& gf, 
					    const VectorInterface& gz) const
  {
    const GlobalVector& f = S.GetGV(gf);
    const GlobalVector& z = S.GetGV(gz);
    double X = f*z;
    
    GlobalVector dz(f.ncomp(),f.n());
    
    dz.zero();

    assert(other);
    assert(dynamic_cast<const AleBaseDiscretization*> (other));    
    PiFSI pi(dynamic_cast<const AleBaseDiscretization*> (other));

    pi.Init(S.GetMesh());
    pi.vmult(dz,z);
    
    ScalarProduct(eta,f,dz);
    
    
    return X;
  }

  double Dwr::ScalarProductWithFluctuationsPrimal(DoubleVector& eta, 
						  const VectorInterface& gf, 
						  const VectorInterface& gz) const
  {
    const GlobalVector& f = S.GetGV(gf);
    const GlobalVector& z = S.GetGV(gz);
    double X = f*z;
    
    GlobalVector dz(f.ncomp(),f.n());
    
    dz.zero();

    assert(other);
    assert(dynamic_cast<const AleBaseDiscretization*> (other));    
    PiFSI pi(dynamic_cast<const AleBaseDiscretization*> (other));

    pi.Init(S.GetMesh());
    pi.vmultprimal(dz,z);
    
    ScalarProduct(eta,f,dz);
    
    return X;
  }
  
  
  /*--------------------------------------------------------*/
  
  DiscretizationInterface* Dwr::CreateOtherDiscretization() const
  {
    DiscretizationInterface* D;
    
    if (S.GetMesh()->dimension()==2) 
      {
	D = new AleDwrQ1Q22d;
      }
    else
      {
	D = new AleDwrQ1Q23d;
      }
    
    D->BasicInit(S.GetParamfile());
    return D;
  }
  
  /*-------------------------------------------------------*/
  
  void Dwr::PrimalResidualsHigher(VectorInterface& gf, const VectorInterface& gu)
  {
    GlobalVector& f = S.GetGV(gf);
    
    f.zero();


    __ESTIMATE = true;
    
    S.Rhs(gf,-0.5);
    S.Form(gf,gu,0.5);
    
    
    
    S.SetDiscretization(*other,true);
    S.Rhs(gf,0.5);
    S.Form(gf,gu,-0.5);

    __ESTIMATE = false;

    S.SetDiscretization(*discretization);
  }
  
  /*--------------------------------------------------------*/
  
  void Dwr::DualResidualsHigher(VectorInterface& gf, 
				const VectorInterface& gu, 
				const VectorInterface& gz, 
				const ProblemDescriptorInterface& PDI)
  {
  }
  
  /*--------------------------------------------------------*/
  
  double Dwr::Estimator(DoubleVector& eta, VectorInterface& gf, 
			VectorInterface& gu, VectorInterface& gz)
  {
    double rho=0, rhostern=0;
    eta.resize(S.GetGV(gz).n());

    other = CreateOtherDiscretization();

    ReInitInterface(other);


    S.SetProblem(*primalproblem);
    PrimalResidualsHigher(gf,gu);
    rho      =  ScalarProductWithFluctuationsPrimal(eta,gf,gz);

    S.SetProblem(*dualproblem);

    other->AddNodeVector("U",&S.GetGV(gu));
    S.AddNodeVector("U",gu);
    __ADJOINT = true;
    PrimalResidualsHigher(gf,gz);
    __ADJOINT = false;
    S.DeleteNodeVector("U");
    other->DeleteNodeVector("U");
    rhostern =  ScalarProductWithFluctuations(eta,gf,gu);




    // Stabilization Part
    S.SetProblem(*primalproblem);
    const Q1Lps2d* LPS2d = dynamic_cast<const Q1Lps2d*> (discretization);
    const Q1Lps3d* LPS3d = dynamic_cast<const Q1Lps3d*> (discretization);
    assert(LPS2d||LPS3d);
    
    GlobalVector& f = S.GetGV(gf);
    f.zero();
    const Equation* EQ = S.GetProblemDescriptor()->GetEquation();
    assert(EQ);
    assert(dynamic_cast<const LpsEquation*> (EQ));

    if (LPS2d) LPS2d->StabForm(f,S.GetGV(gu),*EQ,1.0);
    else if (LPS3d) LPS3d->StabForm(f,S.GetGV(gu),*EQ,1.0);
    else abort();
    

    double RHO_STAB = f * S.GetGV(gz);
    cout << "\t stab:        " << RHO_STAB << endl;
    

    delete other;
    


    
    
    //    std::cout << "P: " << 2.0*rho << "\t D: " << 2.0*rhostern << std::endl;
    
    
    return rho + rhostern + RHO_STAB;
  }
  
  
  /*--------------------------------------------------------*/

  double Dwr::EstimatorEnergy(DoubleVector& eta, VectorInterface& gf, 
				  const VectorInterface& gu)
  {
    eta.resize(S.GetGV(gu).n());
    
    PrimalResidualsHigher(gf,gu);
    return ScalarProductWithFluctuations(eta,gf,gu);
  }
  
/*--------------------------------------------------------*/
}
