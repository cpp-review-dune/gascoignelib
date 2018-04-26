#include "solvers.h"
#include <algorithm>

#include "alediscretization.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "fmatrixblock.h"
#include  "cuthillmckee.h"
#include  "mginterpolatornested.h"

#include "umfilu.h"
#include "sparse_umf.h"

using namespace std;

extern double __DT,__THETA;



namespace Gascoigne
{
  template<int DIM>
  FSISolver<DIM>::FSISolver() : StdSolver(), __IF(0), __IS(0)
  {
  }

  template<int DIM>
  void FSISolver<DIM>::Form(VectorInterface& y, const VectorInterface& x, double d) const
  {
    StdSolver::Form(y,x,d);
  }
  
  
  

  template<int DIM>
  void FSISolver<DIM>::ReInitMatrix() 
  {

   
    if (((_matrixtype=="split") || (_matrixtype=="splitumf"))
	&&(!_directsolver))
      {
	abort();
      }
    else if (((_matrixtype=="fsisplit") || (_matrixtype=="fsisplitumf"))
	     &&(!_directsolver))
      {
	abort();
      }
    
    else StdSolver::ReInitMatrix();    
  }
  
  template<int DIM>
  void FSISolver<DIM>::ComputeIlu(const VectorInterface& gu) const
  {
    if ((_matrixtype=="umfsplit")&&(!_directsolver))
      {
	abort();
      }
    else if (( (_matrixtype=="splitumf")||(_matrixtype=="fsisplitumf") )&&(!_directsolver))
      {
	abort();
      }
    else if ( (_matrixtype=="fsisplit")&&(!_directsolver) )
      {    
	abort();
      }
    else if (!_directsolver)
      {
	int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
	PermutateIlu(gu);
	
	GetIlu()->zero();
	GetIlu()->copy_entries(GetMatrix());
	
	modify_ilu(*GetIlu(),ncomp);
	GetIlu()->compute_ilu();
	
      }
    else StdSolver::ComputeIlu(gu);
  }

  template<int DIM>
  void FSISolver<DIM>::modify_ilu(IluInterface& I,int ncomp) const 
  {      
    if(GetSolverData().GetIluModify().size()==0) return;
    if( GetSolverData().GetIluModify().size()!=ncomp ) {
      cerr << "ERROR: GetSolverData().GetIluModify().size()="<< GetSolverData().GetIluModify().size() << " and ";
      cerr << "ncomp="<< ncomp << endl; 
      abort();
      // assert(GetSolverData().GetIluModify().size()==ncomp);
    }

    for(int c=0;c<ncomp;c++)
      {
	double s = GetSolverData().GetIluModify(c);



	I.modify(c,s);
      }
  }

  
  template<int DIM>
  void FSISolver<DIM>::smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const 
  {
    if ((_matrixtype=="splitumf")&&(!_directsolver))
      {
	abort();
      }
    else if (GetSolverData().GetLinearSmooth()=="ilu")
      {
	double omega = GetSolverData().GetOmega();
	for(int iter=0; iter<niter; iter++)
	  {
	    MatrixResidual(h,x,y);
	    GetIlu()->solve(GetGV(h));
	    Add(x,omega,h);
	  }
      }
    else
      StdSolver::smooth(niter,x,y,h);

  }

  template<int DIM>
  void FSISolver<DIM>::smooth_exact(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const
  {
    StdSolver::smooth_exact(x,y,help);
  }
  


  template<int DIM>
  void FSISolver<DIM>::do_precondition(VectorInterface gx) const
  {
    GlobalVector& H = GetGV(gx);
    assert(H.n() == _precond.size());
    assert(H.ncomp() == _precond[0].size());
    for (int i=0;i<H.n();++i)
      for (int c=0;c<H.ncomp();++c)
	H(i,c) *= _precond[i][c];
  }
  template<int DIM>
  void FSISolver<DIM>::undo_precondition(VectorInterface gx) const
  {
    GlobalVector& H = GetGV(gx);
    assert(H.n() == _precond.size());
    assert(H.ncomp() == _precond[0].size());
    for (int i=0;i<H.n();++i)
      for (int c=0;c<H.ncomp();++c)
	H(i,c) /= _precond[i][c];
  }
  


  /// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  

  template<int DIM>
  MatrixInterface* FSISolver<DIM>::NewMatrix(int ncomp, const std::string& matrixtype)  
  {

#ifdef __WITH_UMFPACK__
    if(_directsolver && _useUMFPACK)             return StdSolver::NewMatrix(ncomp,matrixtype);
#endif
    

    std::string mt = matrixtype;
    cout << mt << endl;
    
    if (matrixtype=="vanka")    mt = "block";
    if (matrixtype=="umfsplit") mt = "block";


    if ((mt=="split") || (mt=="splitumf") )
      {
	abort();
      }
    else if ((mt=="fsisplit") || (mt=="fsisplitumf") )
      {
	abort();
      }
    else 
      return StdSolver::NewMatrix(ncomp,mt);
  }
  
  template<int DIM>	
  IluInterface* FSISolver<DIM>::NewIlu(int ncomp, const string& matrixtype) 
  { 
#ifdef __WITH_UMFPACK__
    //if(_directsolver && _useUMFPACK)             return new UmfIluLong(GetMatrix());
    if(_directsolver && _useUMFPACK)             return new UmfIlu(GetMatrix());
    //  if(_directsolver && _useUMFPACK) {cout<<"keine ahnung was hier passiert"<<endl; abort();}
#endif
    cout << matrixtype << endl;
    
    if (matrixtype=="split")
      {
	abort();
      }
    else if (matrixtype=="splitumf")
      {
	abort();
      }
    else if (matrixtype=="fsisplitumf")
      {
	abort();
      }
    else if (matrixtype=="fsisplit")
      {
	abort();
      }
    else return StdSolver::NewIlu(ncomp, matrixtype);
    
    // SplittingIlu<DIM>* TMP = new SplittingIlu<DIM>;
    // TMP->SetInterface(GetAleDiscretization()->GetFluidL2G(), 
    // 		      GetAleDiscretization()->GetSolidL2G(), 
    // 		      GetAleDiscretization()->GetFluidG2L(), 
    // 		      GetAleDiscretization()->GetSolidG2L(),
    // 		      GetAleDiscretization()->GetInterfaceNodes());
    
    // return TMP;
  }
  
  
  
  /*-------------------------------------------------------*/
  
  template<int DIM>
  void FSISolver<DIM>::DeleteSolidPressure(VectorInterface& gf) const
  {
    ////////////////////
    
    const HASHSET<int>&     i_nodes = GetAleDiscretization()->GetInterfaceNodes();
    const vector<int>&      s_nodes = GetAleDiscretization()->GetSolidL2G();
    const vector<int>&      f_nodes = GetAleDiscretization()->GetFluidL2G();
    /*
      for (auto node : s_nodes)
      if (i_nodes.find(node)==i_nodes.end())
      {
      GetGV(gf)(node,0) = 0.0;
      }
    */
    if(GetProblemDescriptor()->GetName()=="solid_euler")
	{
		//Delete Pressure	
		for (vector<int>::const_iterator it = s_nodes.begin();it!=s_nodes.end();++it)
		{
		 	 GetGV(gf)(*it,0) = 0.0;
		 	
		}
		for (vector<int>::const_iterator it = f_nodes.begin();it!=f_nodes.end();++it)
		{
		 if (i_nodes.find(*it)==i_nodes.end())
		  	{
		  		GetGV(gf)(*it,0) = 0.0;
		  		//GetGV(gf)(*it,1) = 0.0;
		  	    //GetGV(gf)(*it,2) = 0.0;
		  		//GetGV(gf)(*it,3) = 0.0;
		  	}
		  	
		}
	}

	
	if(GetProblemDescriptor()->GetName()=="fluid_stat")
	{
		//Delete Pressure	
		for (vector<int>::const_iterator it = s_nodes.begin();it!=s_nodes.end();++it)
		  if (i_nodes.find(*it)==i_nodes.end())
		{
		  GetGV(gf)(*it,0) = 0.0;
		}
	}
	//Delete Components in f on solid domain
	if(GetProblemDescriptor()->GetName()=="fluid_stat")
	{
		for (vector<int>::const_iterator it = s_nodes.begin();it!=s_nodes.end();++it)
		{
		
		  	GetGV(gf)(*it,1) = 0.0;
		  	GetGV(gf)(*it,2) = 0.0;
		  	GetGV(gf)(*it,3) = 0.0;
		}
	}
	
	if(GetProblemDescriptor()->GetName()=="solid_euler")
	{
	// PointVisu("residuum",GetGV(gf),0);

	}
    // for (auto node : s_nodes)
    //   for (int c=0;c<DIM;++c)
    // 	GetGV(gf)(node,c+1+DIM) = 0.0;
    if(GetProblemDescriptor()->GetName()!="fluid_stat" ||GetProblemDescriptor()->GetName()!="solid_euler")
	{
		if(GetSolverLabel()=="fsi_reduced")
		{
			//Delete Pressure	
			for (vector<int>::const_iterator it = s_nodes.begin();it!=s_nodes.end();++it)
			  if (i_nodes.find(*it)==i_nodes.end())
			{
			  GetGV(gf)(*it,0) = 0.0;
			} 
		}
		if(GetSolverLabel()=="def_solid")
		{

			for (vector<int>::const_iterator it = f_nodes.begin();it!=f_nodes.end();++it)
			  if (i_nodes.find(*it)==i_nodes.end())
				{
				  for(int i=0;i<DIM;i++)	
				  	GetGV(gf)(*it,i) = 0.0;
				} 
		}
		if(GetSolverLabel()=="meshmotion")
		{
			for (vector<int>::const_iterator it = s_nodes.begin();it!=s_nodes.end();++it)
			{
				for(int i=0;i<DIM;i++)	
				  	GetGV(gf)(*it,i) = 0.0;
			}
		}
	}
    
  }

  template<int DIM>
  void FSISolver<DIM>::SetBoundaryVectorZero(VectorInterface& gf) const
  {
    StdSolver::SetBoundaryVectorZero(gf);
  
    DeleteSolidPressure(gf);
  }

  template<int DIM>
  void FSISolver<DIM>::SetBoundaryVector(VectorInterface& gf) const
  {
    StdSolver::SetBoundaryVector(gf);
  
    DeleteSolidPressure(gf);
  }


  template<int DIM>
  void FSISolver<DIM>::AssembleMatrix(const VectorInterface& gu, double d)
  {
    
    StdSolver::AssembleMatrix(gu,d);
    

    if ((_directsolver)||(_matrixtype=="block"))
      {

	// Modify for pressure zero in Solid-Part
	const HASHSET<int>&     i_nodes = GetAleDiscretization()->GetInterfaceNodes();
	const vector<int>&      s_nodes = GetAleDiscretization()->GetSolidL2G();
	const vector<int>&      f_nodes = GetAleDiscretization()->GetFluidL2G();
	
	if(GetProblemDescriptor()->GetName()=="fluid_stat")
	{
		vector<int> cv1,cv2;
		cv1.push_back(0);
		cv1.push_back(1);
		cv1.push_back(2);
		cv1.push_back(3);

		cv2.push_back(1);
		cv2.push_back(2);
		cv2.push_back(3);
		for (int i=0;i<s_nodes.size();++i)
		  {
		  	if (i_nodes.find(s_nodes[i])==i_nodes.end())
		   		 GetMatrix()->dirichlet(s_nodes[i],cv1);
		   	else
		   		 GetMatrix()->dirichlet(s_nodes[i],cv2);
		  }
	}
	else if(GetProblemDescriptor()->GetName()=="solid_euler")
	{
		vector<int> cv1;
		cv1.push_back(0);
		cv1.push_back(1);
		cv1.push_back(2);
		cv1.push_back(3);
		
	    vector<int> cv2;
	    cv2.push_back(0);
			
		for (int i=0;i<f_nodes.size();++i)
		  {
		  	if (i_nodes.find(f_nodes[i])==i_nodes.end())
		   		 GetMatrix()->dirichlet(f_nodes[i],cv1);
		  }
		for (int i=0;i<s_nodes.size();++i)
		 {
		    	GetMatrix()->dirichlet(s_nodes[i],cv2);
		 } 
	}
	else
	{
	//Delete Pressure variable in Matrix
		if(GetSolverLabel()=="fsi_reduced"||GetSolverLabel()=="fsi_main")
		{
			vector<int> cv;
			cv.push_back(0);
			for (int i=0;i<s_nodes.size();++i)
		  		if (i_nodes.find(s_nodes[i])==i_nodes.end())
		    	GetMatrix()->dirichlet(s_nodes[i],cv);
		}
		if(GetSolverLabel()=="def_solid")
		{

			vector<int> cv;
			for(int i=0;i<DIM;i++) cv.push_back(i);

			for (int i=0;i<f_nodes.size();++i)
		  	{	
		  		if (i_nodes.find(f_nodes[i])==i_nodes.end())
		    		GetMatrix()->dirichlet(f_nodes[i],cv);
		    }
		}
		if(GetSolverLabel()=="meshmotion")
		{
			vector<int> cv;
			for(int i=0;i<DIM;i++) cv.push_back(i);
	
			for (int i=0;i<s_nodes.size();++i)
		  		//if (i_nodes.find(s_nodes[i])==i_nodes.end())
		    	GetMatrix()->dirichlet(s_nodes[i],cv);
		}
	}
	// cv.clear();
	// for (int i=0;i<DIM;++i) cv.push_back(i+1+DIM);
	// for (int node=0;node<GetMesh()->nnodes();++node)
	//   {
	//       GetMatrix()->dirichlet_only_row(node,cv);
	//   }
      }
  }


  template<int DIM>
  DiscretizationInterface* FSISolver<DIM>::NewDiscretization(int dimension, const string& discname)
  {
    if (dimension==2)
      {
	if      (discname=="AleQ1")             {abort();  return new AleQ12d;}
	//	else if (discname=="AleQ2")               return new AleQ22d;
	else if (discname=="AleQ1Lps")          {abort(); return new AleQ1Lps2d;}
	
	else if (discname=="AleQ2Lps")
	  {
	    DiscretizationInterface* X = new AleQ2Lps2d;
	    //vector<int> delf,dels;
	    //	    for (int i=0;i<DIM;++i) delf.push_back(i+1+DIM);
	    //dels.push_back(0);
	    //dynamic_cast<AleQ2Lps2d*>(X)->InitInterfaceComponents(delf,dels);
	    return X;
	  }
	
	else return StdSolver::NewDiscretization(dimension, discname);
      }
    else if (dimension==3)
      {
	// if      (discname=="AleQ1")               return new AleQ13d;
	if (discname=="AleQ1Lps")            return new AleQ1Lps3d;
	else if (discname=="AleQ2Lps")            
	  {
	    DiscretizationInterface* X = new AleQ2Lps3d;
	    //vector<int> delf,dels;
	    //	    for (int i=0;i<DIM;++i) delf.push_back(i+1+DIM);
	    //dels.push_back(0);
	    //dynamic_cast<AleQ2Lps3d*>(X)->InitInterfaceComponents(delf,dels);
	    return X;
	  }
	
	// else if (discname=="AleQ2Lps")            return new AleQ2Lps3d;
	// else
	return StdSolver::NewDiscretization(dimension, discname);
      }
    else abort();
  }


  template<int DIM>
  void FSISolver<DIM>::reinit_interface_element(int en, const nvector<int>& indices, 
						HASHMAP<int, std::vector<int> >& solid_interface_cells, 
						HASHMAP<int, std::vector<int> >& fluid_interface_cells,
						HASHSET<int> & interface_nodes, int material)
  {
    vector<int> ni;
    for (int i=0;i<indices.size();++i)
      {
	if(interface_nodes.find(indices[i])!=interface_nodes.end())
	  {
	    ni.push_back(i);
	  }
      }
  
    if (ni.size()>0)
      {
	if(material==1) 
	  {//solid cell     
	    solid_interface_cells[en]=ni;
	  }
	else if(material==2) 
	  {//fluid cell
	    fluid_interface_cells[en]=ni;
	  }
	else
	  {cout<<"	Fluid Cells need to have the material value 2 and solid cells the material value 1!" <<endl; abort();}
			
      }
  }
     			    
  template<int DIM>
  void FSISolver<DIM>::reinit_element(int en, const nvector<int>& indices, 
				      HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
				      set<int>& fluid_nodes, set<int>& solid_nodes,int material)
  {
    
    if(material==1) 
      {//solid cell
	for (int i=0;i<indices.size();++i)
	  { 
	    solid_nodes.insert(indices[i]);
	    solid_cells.insert(en);
	  }
      }
    else if(material==2) 
      {//fluid cell
	for (int i=0;i<indices.size();++i)
	  { 
	    fluid_nodes.insert(indices[i]);
	    fluid_cells.insert(en);
	  }
      }
    else
      {cout<<"	Fluid Cells need to have the material value 2 and solid cells the material value 1!" <<endl; abort();}
		

  }


  template<int DIM>
  void FSISolver<DIM>::ReInitInterface(AleBaseDiscretization* ALEDISC)
  {
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
    
    
    int dim = GetMesh()->dimension();

    if (dim==2)
      {	
	const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*> (GetMesh());
	assert(M);
	if ((GetDiscretization()->GetName()=="Q1 Ale 2d Lps")||
	    (GetDiscretization()->GetName()=="Q1 Ale 2d"))
	  {
	    // Einsortieren der Solid und Fluid Nodes
	    for (int c=0;c<M->ncells();++c)
	      {
		reinit_element(c, M->IndicesOfCell(c), fluid_cells, solid_cells, fluid_nodes, solid_nodes,M->material(c));
	      }
	    // Interface Node: Sowohl Fluid als auch Solid Node
	    // Kann erst aufgerufen werden wenn man durch alle Zellen einmal durch ist!
	    for(set<int>::const_iterator it = fluid_nodes.begin();it!=fluid_nodes.end();++it)  
	      if(solid_nodes.find(*it)!=solid_nodes.end()) interface_nodes.insert(*it);
	    // Interface Cells und Interfacenodes on InterfaceCells abspeichern	
	    for (int c=0;c<M->ncells();++c)
	      {		
		reinit_interface_element(c, M->IndicesOfCell(c), solid_interface_cells, 
					 fluid_interface_cells, interface_nodes,M->material(c) );
	      }
	  }
	else if ((GetDiscretization()->GetName()=="Q2 Ale 2d Lps")||
		 (GetDiscretization()->GetName()=="Q2 Ale 2d"))
	  {
	    // Einsortieren der Solid und Fluid Nodes
	    for (int c=0;c<M->npatches();++c)
	      {
		reinit_element(c, *(M->IndicesOfPatch(c)), fluid_cells, solid_cells, fluid_nodes, solid_nodes,M->material_patch(c));
	      }
	    // Interface Node: Sowohl Fluid als auch Solid Node
	    // Kann erst aufgerufen werden wenn man durch alle Zellen einmal durch ist!
	    for(set<int>::const_iterator it = fluid_nodes.begin();it!=fluid_nodes.end();++it)  
	      if(solid_nodes.find(*it)!=solid_nodes.end()) interface_nodes.insert(*it);
	    // Interface Cells und Interfacenodes on InterfaceCells abspeichern	
	    for (int c=0;c<M->npatches();++c)
	      {		
		reinit_interface_element(c, *(M->IndicesOfPatch(c)), solid_interface_cells, 
					 fluid_interface_cells, interface_nodes,M->material_patch(c) );
	      }
	  }
	else abort();
      
      }
    else if (dim==3)
      {	
	const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (GetMesh());
	assert(M);
	
	if (GetDiscretization()->GetName()=="Q1 Ale 3d Lps")
	  {
	    // Einsortieren der Solid und Fluid Nodes
	    for (int c=0;c<M->ncells();++c)
	      {
		reinit_element(c, M->IndicesOfCell(c), fluid_cells, solid_cells, fluid_nodes, solid_nodes,M->material(c));
	      }
	    // Interface Node: Sowohl Fluid als auch Solid Node
	    // Kann erst aufgerufen werden wenn man durch alle Zellen einmal durch ist!
	    for(set<int>::const_iterator it = fluid_nodes.begin();it!=fluid_nodes.end();++it)  
	      if(solid_nodes.find(*it)!=solid_nodes.end()) interface_nodes.insert(*it);
	    // Interface Cells und Interfacenodes on InterfaceCells abspeichern	
	    for (int c=0;c<M->ncells();++c)
	      {		
		reinit_interface_element(c, M->IndicesOfCell(c), solid_interface_cells, 
					 fluid_interface_cells, interface_nodes,M->material(c) );
	      }
	  }
	else if (GetDiscretization()->GetName()=="Q2 Ale 3d Lps")
	  {
	    // Einsortieren der Solid und Fluid Nodes
	    for (int c=0;c<M->npatches();++c)
	      {
		reinit_element(c, *(M->IndicesOfPatch(c)), fluid_cells, solid_cells, fluid_nodes, solid_nodes,M->material_patch(c));
	      }
	    // Interface Node: Sowohl Fluid als auch Solid Node
	    // Kann erst aufgerufen werden wenn man durch alle Zellen einmal durch ist!
	    for(set<int>::const_iterator it = fluid_nodes.begin();it!=fluid_nodes.end();++it)  
	      if(solid_nodes.find(*it)!=solid_nodes.end()) interface_nodes.insert(*it);
	    // Interface Cells und Interfacenodes on InterfaceCells abspeichern	
	    for (int c=0;c<M->npatches();++c)
	      {		
		reinit_interface_element(c, *(M->IndicesOfPatch(c)), solid_interface_cells, 
					 fluid_interface_cells, interface_nodes,M->material_patch(c) );
	      }
	  }
	else 
	  {
	    std::cout << GetDiscretization()->GetName() << std::endl;
					
	    abort();
	  }
	

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

  // --------------------------------------------------

  template<int DIM>
  void FSISolver<DIM>::NewMesh(const MeshInterface* mp)
  {
    StdSolver::NewMesh(mp);

    ReInitInterface(GetAleDiscretization());
  }

  
  template<>
  void FSISolver<2>::PointVisu(const string& name, const GlobalVector& u, int iter) const
  {
	cout<<	"PointVisu for DIM=2 not written "<<endl; abort();
    StdSolver::PointVisu(name,u,iter);
  }
  template<>
  void FSISolver<3>::PointVisu(const string& name, const GlobalVector& u, int iter) const
  {

		/*VectorInterface def("def");
		const GlobalVector& DEF = GetGV(def);
		
		GlobalVector U;
		if(GetProblemDescriptor()->GetName()=="fluid_stat")
    	{
    	 U.ncomp()  = u.ncomp()+1;
		 //1 pressure ,3 fluid ,1 domain
		 assert(1+3+1==U.ncomp());
    	}
    	else if(GetProblemDescriptor()->GetName()=="solid_euler")
    	{
    		U.ncomp()  = u.ncomp()+1;
			//1?,3 solid ,1 domain
			assert(3+1+1==U.ncomp());
    	}
    	else 
    	{
    	 U.ncomp()  = 2*u.ncomp();
		 //1 pressure ,3 fluid ,3 solid ,1 domain
		 assert(2*3+2==U.ncomp());
    	}
		U.resize(u.n());
		
		const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (GetMesh());
		assert(M);
		
		for (int i=0;i<u.n();++i)
		  {
			const Vertex3d v = M->vertex3d(i);
			//int domain = chi(v);
			int domain;
			if (find (GetAleDiscretization()->GetFluidL2G().begin(), GetAleDiscretization()->GetFluidL2G().end(), i) != GetAleDiscretization()->GetFluidL2G().end())
			  {
				if (GetAleDiscretization()->GetInterfaceNodes().find(i)!= GetAleDiscretization()->GetInterfaceNodes().end())
				  domain=0;
				else
				  domain=-1;
			  }
			else
			  {
				if (find (GetAleDiscretization()->GetSolidL2G().begin(), GetAleDiscretization()->GetSolidL2G().end(), i) != GetAleDiscretization()->GetSolidL2G().end())
				  {
				if (GetAleDiscretization()->GetInterfaceNodes().find(i)!= GetAleDiscretization()->GetInterfaceNodes().end())
				  domain=0;
				else
				  domain=1;
				  }
				else {cout<<"error writing mesh. Node neither fluid nor solid"<<endl; abort();}
			  }
	
			for (int c=0;c<u.ncomp();++c)
			  U(i,c) = u(i,c);
			if(GetProblemDescriptor()->GetName()=="fsi")  
			{
				for (int c=0;c<3;++c)
			  		U(i,c+3+1) = DEF(i,c);
			}
			U(i,U.ncomp()-1) = domain;
		  }
		  		StdSolver::PointVisu(name,U,iter);
		  */
		StdSolver::PointVisu(name,u,iter);
		
		
	   }
  //////////////////////////////////////////////////


  
	  

  

  template class FSISolver<2>;
  template class FSISolver<3>;
}
