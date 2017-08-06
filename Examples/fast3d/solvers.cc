#include "solvers.h"

#include "alediscretization.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
//#include "splittingilu.h"
#include "umfilu.h"
#include "fmatrixblock.h"
#include  "cuthillmckee.h"
#include  "mginterpolatornested.h"

using namespace std;


int COMPILU;


namespace Gascoigne
{
  template<int DIM>
  FSISolver<DIM>::FSISolver() : StdSolver(), __IF(0), __IS(0)
  {
    COMPILU = 0;
    cout << "CI: " << COMPILU << endl;

  }

  template<int DIM>
  void FSISolver<DIM>::Form(VectorInterface& y, const VectorInterface& x, double d) const
  {
    StdSolver::Form(y,x,d);
  }
  
  
  /// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ILU-STUFF

  template<int DIM>
  void FSISolver<DIM>::ReInitMatrix() 
  {
    COMPILU = 0;
    cout << "CI: " << COMPILU << endl;
    
    if ((_matrixtype=="umfsplit")&&(!_directsolver))
      {
	abort();
	
	GetDiscretization()->InitFilter(GetPfilter());
	SparseStructure SA;
	GetDiscretization()->Structure(&SA);
	
	if (GetFaceDiscretization())
	  GetFaceDiscretization()->Structure(&SA);
	
	AddPeriodicNodes(&SA);
	
	GetMatrix()->ReInit(&SA);
	
	assert(__IF); __IF->ReInit(&SA);
	assert(__IS); __IS->ReInit(&SA);
      }
    else if (((_matrixtype=="split") || (_matrixtype=="splitumf"))
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
  void FSIMultiLevelSolver<DIM>::ComputeIlu(VectorInterface& u)
  {
    //    COMPILU++;
    //    cout << "CI: " << COMPILU << endl;//
    //    if ((COMPILU%4)==1) 
StdMultiLevelSolver::ComputeIlu(u);
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
  void FSISolver<DIM>::Precondition(const IluInterface& M,GlobalVector &x) const
  {
    abort();
  }
  

  template<int DIM>
  void FSISolver<DIM>::Precondition_SINGLE(const IluInterface& M,GlobalVector &x) const
  {
    abort();
  }

  template<int DIM>
  int  FSISolver<DIM>::BiCGSTAB(const MatrixInterface &A, 
				GlobalVector &x, const GlobalVector &b,
				const IluInterface &M, 
				int &max_iter, double &tol) const 
  {
    abort();
  }



  template<int DIM>
  int  FSISolver<DIM>::BiCGSTAB_SINGLE(const MatrixInterface &A, 
				       GlobalVector &x, const GlobalVector &b,
				       const IluInterface &M, 
				       int &max_iter, double &tol) const 
  {
    int nc = x.ncomp();
    int nn = x.n();
    
    double resid=0;
    double rho_1=0, rho_2=0, alpha=0, beta=0, omega=0;
    GlobalVector p(nc,nn), phat(nc,nn), s(nc,nn), shat(nc,nn), t(nc,nn), v(nc,nn);
    
    double normb = b.norm();
    
    GlobalVector r = b;
    A.vmult(r,x,-1.0);
    GlobalVector rtilde = r;
    
    if (normb == 0.0)
      normb = 1;
    
    resid = r.norm()/normb;
    tol *= resid;
    
    if (resid <= tol) 
      {
	tol = resid;
	max_iter = 0;
	return 0;
      }
    
    for (int i = 1; i <= max_iter; i++) 
      {
	rho_1 = rtilde * r;
	if (rho_1 == 0) 
	  {
	    tol = r.norm() / normb;
	    std::cerr << "BiCgStab exploded, <rtilde,r>=0" << std::endl;
	    return 2;
	  }

	p = r;
	if (i > 1)
	  {
	    beta = (rho_1/rho_2) * (alpha/omega);
	    p.add(beta,p, -omega, v);
	  }
	
	phat = p;
	
	///////////////////// Vorkonditionierung!!!!!!
	//	Precondition_SINGLE(M,phat);
	/////////////////////
	
	v.zero();
	A.vmult(v,phat,1.0);
	alpha = rho_1 / (rtilde * v);
	s=r;
	s.add(-alpha,v);
	resid = s.norm()/normb;

	//	cout << "\t -> BiCGStab\t " << i << "\t" << resid << endl;
	cout << i << " " << resid << endl;
	if (resid< tol) 
	  {
	    x.add(alpha,phat);
	    tol = resid;
	    return 0;
	  }
	shat = s;

	/////////////////// Vorkonditionieren

	//	Precondition_SINGLE(M,shat);

	/////////////////// Vorkonditionieren

	t.zero();
	A.vmult(t,shat,1.0);
	
	omega = (t*s) / (t*t);
	x.add(alpha,phat,omega,shat);
	r = s;
	r.add (-omega,t);
	
	rho_2 = rho_1;
	resid = r.norm() / normb;
	cout << i << " " << resid << endl;
	//	cout << "\t -> BiCGStab\t " << i << "\t" << resid << endl;
	if (resid < tol) 
	  {
	    tol = resid;
	    max_iter = i;
	    return 0;
	  }
	if (omega == 0) 
	  {
	    tol = r.norm()/normb;
	    std::cerr << "BiCgStab STOP, <rtilde,r>=0" << std::endl;
	    return 3;
	  }
      }
    
    tol = resid;
    return 1;
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
    if(_directsolver && _useUMFPACK)             return new UmfIluLong(GetMatrix());
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
    vector<int> cv;
    for (int i=0;i<s_nodes.size();++i)
      if (i_nodes.find(s_nodes[i])==i_nodes.end())
	{
	  GetGV(gf)(s_nodes[i],0) = 0.0;
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
    
    
    if ((_directsolver)||(_matrixtype=="block")||(_matrixtype=="sparse"))
      {

	// Modify for pressure zero in Solid-Part
	const HASHSET<int>&     i_nodes = GetAleDiscretization()->GetInterfaceNodes();
	const vector<int>&      s_nodes = GetAleDiscretization()->GetSolidL2G();
	vector<int> cv;
	cv.push_back(0);
	for (int i=0;i<s_nodes.size();++i)
	  if (i_nodes.find(s_nodes[i])==i_nodes.end())
	    GetMatrix()->dirichlet(s_nodes[i],cv);
      }


    // cout << _matrixtype << endl;
  }


  template<int DIM>
  DiscretizationInterface* FSISolver<DIM>::NewDiscretization(int dimension, const string& discname)
  {
    if (dimension==2)
      {
	if      (discname=="AleQ1")               return new AleQ12d;
	//	else if (discname=="AleQ2")               return new AleQ22d;
	else if (discname=="AleQ1Lps")            return new AleQ1Lps2d;
	else if (discname=="AleQ2Lps")            return new AleQ2Lps2d;
	else return StdSolver::NewDiscretization(dimension, discname);
      }
    else if (dimension==3)
      {
	// if      (discname=="AleQ1")               return new AleQ13d;
	if (discname=="AleQ1Lps")            return new AleQ1Lps3d;
	else if (discname=="AleQ2Lps")            return new AleQ2Lps3d;
	// else if (discname=="AleQ2Lps")            return new AleQ2Lps3d;
	// else
	return StdSolver::NewDiscretization(dimension, discname);
      }
    else abort();
  }

  template<>
  void FSISolver<2>::reinit_element(int en, const nvector<int>& indices, 
				    HASHMAP<int, std::vector<int> >& solid_interface_cells, 
				    HASHMAP<int, std::vector<int> >& fluid_interface_cells,
				    HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
				    HASHSET<int> & interface_nodes,
				    set<int>& fluid_nodes, set<int>& solid_nodes)
  {
    Chi chi;
    
    const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*> (GetMesh());
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
	cerr.precision(20);
	
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

  template<>
  void FSISolver<3>::reinit_element(int en, const nvector<int>& indices, 
				    HASHMAP<int, std::vector<int> >& solid_interface_cells, 
				    HASHMAP<int, std::vector<int> >& fluid_interface_cells,
				    HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
				    HASHSET<int> & interface_nodes,
				    set<int>& fluid_nodes, set<int>& solid_nodes)
  {
    Chi chi;
    
    
    const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (GetMesh());
    assert(M);

    int nf=0;
    int ns=0;
    
    vector<int> ni;
    //    cout << endl << endl;
    
    for (int i=0;i<indices.size();++i)
      {
	int domain = chi(M->vertex3d(indices[i]));
	//	cout << domain << "\t" << M->vertex3d(indices[i]) << endl;
	
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
	  for (int c=0;c<M->ncells();++c)
	    reinit_element(c, M->IndicesOfCell(c), solid_interface_cells, 
			   fluid_interface_cells, fluid_cells, solid_cells, 
			   interface_nodes, fluid_nodes, solid_nodes);
      
	else if ((GetDiscretization()->GetName()=="Q2 Ale 2d Lps")||
		 (GetDiscretization()->GetName()=="Q2 Ale 2d"))
	  for (int c=0;c<M->npatches();++c)
	    reinit_element(c, *(M->IndicesOfPatch(c)), solid_interface_cells, 
			   fluid_interface_cells, fluid_cells, solid_cells, 
			   interface_nodes, fluid_nodes, solid_nodes);
	else abort();
      
      }
    else if (dim==3)
      {	
	const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (GetMesh());
	assert(M);
	
	if (GetDiscretization()->GetName()=="Q1 Ale 3d Lps")
	  for (int c=0;c<M->ncells();++c)
	    reinit_element(c, M->IndicesOfCell(c), solid_interface_cells, fluid_interface_cells, fluid_cells, solid_cells, 
			   interface_nodes, fluid_nodes, solid_nodes);
	else if (GetDiscretization()->GetName()=="Q2 Ale 3d Lps")
	  for (int c=0;c<M->npatches();++c)
	    reinit_element(c, *(M->IndicesOfPatch(c)), solid_interface_cells, fluid_interface_cells, fluid_cells, solid_cells, 
			   interface_nodes, fluid_nodes, solid_nodes);
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

    GlobalVector U;
    U.ncomp()  = u.ncomp()+1;
    assert(2*2+2==U.ncomp());
    U.resize(u.n());
    
    
    const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*> (GetMesh());
    assert(M);
    Chi chi;
    
    for (int i=0;i<u.n();++i)
      {
	const Vertex2d v = M->vertex2d(i);
	int domain = chi(v);
	for (int c=0;c<u.ncomp();++c)
	  U(i,c) = u(i,c);
	U(i,u.ncomp()) = domain;
      }
    StdSolver::PointVisu(name,U,iter);
  }
  template<>
  void FSISolver<3>::PointVisu(const string& name, const GlobalVector& u, int iter) const
  {
    GetDiscretization()->HNAverage(const_cast<GlobalVector&>(u)); 
    GlobalVector U;
    U.ncomp()  = u.ncomp()+1;
    assert(2*3+2==U.ncomp());
    U.resize(u.n());
    
    const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (GetMesh());
    assert(M);
    Chi chi;
    
    for (int i=0;i<u.n();++i)
      {
	const Vertex3d v = M->vertex3d(i);
	int domain = chi(v);
	for (int c=0;c<u.ncomp();++c)
	  U(i,c) = u(i,c);
	U(i,u.ncomp()) = domain;
      }
    StdSolver::PointVisu(name,U,iter);
  }
  



  template class FSISolver<2>;
  template class FSISolver<3>;
  template class FSIMultiLevelSolver<2>;
  template class FSIMultiLevelSolver<3>;
}
