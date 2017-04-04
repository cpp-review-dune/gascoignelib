#include "alesolver.h"
#include "aleinterpolator.h"
#include "alediscretization.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "umfilu.h"
#include "sparse_umf.h"
#include "sparseblockmatrix.h"
#include "fmatrixblock.h"
#include  "cuthillmckee.h"
#include  "sparseblockilu.h"
#include "aleblocks.h"


using namespace std;


namespace Gascoigne
{
  
  StopWatch S1,S2,S25,S3,S4,S5,S6,S_LIN;
  

  

   int  AleSolver::Richardson(const MatrixInterface &A, 
			      GlobalVector &x, const GlobalVector &b,
			      const IluInterface &M, 
			      int &max_iter, 
			      double &tol) const 
   {
     GlobalVector H(x.ncomp(), x.n());
     GlobalVector R(x.ncomp(), x.n());
     GlobalVector HH(x.ncomp(), x.n());
     
     x.zero();
     double start_res=0;
     double res      =0;
     int i;
     for (i=0;i<max_iter;++i)
       {
	 R=b;
	 if (i>0)
	   A.vmult(R,x,-1.0);
	 else
	   start_res = R.norm();
	 
	 res = R.norm();

	 if (fabs(res)<1.e-18)    break;
	 if (res/start_res < tol) break;
	 
	 
	 
	 H=R;
	 M.solve(H);
	 
	 HH.zero();
	 A.vmult(HH,H,1.0);
	 double s = (HH*R)/ (HH*HH);
	 x.add(s,H);
	 
       }
     cout << "R: " << i << "\t" << start_res << "\t " << res << "\t"
	  << res/start_res << endl;

     return 0;
     
   }
  
  
  int  AleSolver::BiCGSTAB(const MatrixInterface &A, 
			   GlobalVector &x, const GlobalVector &b,
			   const IluInterface &M, 
			   int &max_iter, 
			   double &tol) const 
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
    
    //    cout << "\t -> BiCGStab\t " << 0 << "\t" << resid << endl;
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
	    return 2;
	  }
	if (i == 1)
	  p = r;
	else 
	  {
	    beta = (rho_1/rho_2) * (alpha/omega);
	    p=r;
	    p.add(beta,p, -omega, v);
	  }
	phat = p;
	M.solve(phat);
	v.zero();
	A.vmult(v,phat,1.0);
	alpha = rho_1 / (rtilde * v);
	s=r;
	s.add(-alpha,v);
	resid = s.norm()/normb;
	//	cout << "\t -> BiCGStab\t " << i << "\t" << resid << endl;
	if (resid< tol) 
	  {
	    x.add(alpha,phat);
	    tol = resid;
	    return 0;
	  }
	shat = s;
	M.solve(shat);
	t.zero();
	A.vmult(t,shat,1.0);
	
	omega = (t*s) / (t*t);
	x.add(alpha,phat,omega,shat);
	r = s;
	r.add (-omega,t);
	
	rho_2 = rho_1;
	resid = r.norm() / normb;
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
	    return 3;
	  }
      }
    
    tol = resid;
    return 1;
  }







  

  AleSolver::AleSolver() : StdSolver(), 
			   __DUALMATRIX(false),
			   __splittingsmoother(false),
			   __splitting_fluid_exact(true),
			   __splitting_solid_exact(true),
			   __AF(0), __AS(0), __IF(0), __IS(0)
  {
  }  

  // ==================================================
  
  AleSolver::~AleSolver()
  {
    if (__AF) delete __AF; __AF = 0;
    if (__AS) delete __AS; __AS = 0;
    if (__IF) delete __IF; __IF = 0;
    if (__IS) delete __IS; __IS = 0;
  }
  
  // ==================================================

  void AleSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
  {
    S3.start();
    StdSolver::Form(gy,gx,d);

    S3.stop();
  }

  
  template<class FF, class FS, class FA>
  void AleSolver::SplitMatrix(SparseBlockMatrix<FF>* AF,SparseBlockMatrix<FS>* AS, const SparseBlockMatrix<FA>* A)
  {
    assert(A);
    assert(AS);
    assert(AF);

    const vector<int>& fluid_l2g = GetAleDiscretization()->GetFluidL2G();
    const vector<int>& solid_l2g = GetAleDiscretization()->GetSolidL2G();
    const HASHMAP<int,int>& fluid_g2l = GetAleDiscretization()->GetFluidG2L();
    const HASHMAP<int,int>& solid_g2l = GetAleDiscretization()->GetSolidG2L();
    const HASHSET<int>& interface_nodes = GetAleDiscretization()->GetInterfaceNodes();
	
    S6.start();
    
    if (solid_l2g.size()>0)
      CopySubMatrixSolid(AS,A, solid_l2g, solid_g2l,interface_nodes);
    
    assert(fluid_l2g.size()>0);
    CopySubMatrixFluid(AF,A, fluid_l2g, fluid_g2l, interface_nodes);

    S6.stop();
	
    
    S2.stop();
    
    // Compute-Ilu's
    S5.start();
    int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
    assert(__IS);
    assert(__IF);
    
    if (solid_l2g.size()>0)
      {
    	if (!__splitting_solid_exact)
    	  {
    	    nvector<int> PS(solid_l2g.size());
    	    if (1)
    	      {
    		CuthillMcKee cmc(__AS->GetStencil());
    		cmc.Permutate(PS);
    	      }
    	    __IS->ConstructStructure(PS,*__AS);
    	    __IS->zero();
    	    __IS->copy_entries(__AS);
    	    modify_ilu(*__IS,ncomp);
    	    __IS->compute_ilu();
    	  }
    	else
    	  {
    	    ////////////////////////////// WITH SparseUMF
    	    vector<int> P;
    	    __IS->ConstructStructure(P,*__AS);
    	    __IS->zero();
    	    __IS->copy_entries(__AS);
    	    __IS->compute_ilu();
    	  }
      }
    
    
    
    if (!__splitting_fluid_exact)
      {
    	nvector<int> PF(fluid_l2g.size());
    	if (1)
    	  {
    	    CuthillMcKee cmc(__AF->GetStencil());
    	    cmc.Permutate(PF);
    	  }
    	__IF->ConstructStructure(PF,*__AF);
    	__IF->zero();
    	__IF->copy_entries(__AF); 

	for(int c=0;c<ncomp;c++)
	  {
	    double s = 0.0;
	    if ((c>=1)&&(c<=3)) s = 0.1;
	    if (c==0) s = 0.01;
	    __IF->modify(c,s);
	  }
	//    	modify_ilu(*__IF,ncomp);

    	__IF->compute_ilu();
    	S5.stop();
      }
    else
      {
    	////////////////////////////// WITH SparseUMF
    	vector<int> P;
    	__IF->ConstructStructure(P,*__AF);
    	__IF->zero();
    	__IF->copy_entries(__AF);
    	__IF->compute_ilu();
      }
    
  }
  
/*-------------------------------------------------------*/

  void AleSolver::SetBoundaryVectorZero(VectorInterface& gf) const
  {
    StdSolver::SetBoundaryVectorZero(gf);
    
    const HASHSET<int>&     i_nodes = GetAleDiscretization()->GetInterfaceNodes();
    const vector<int>&      s_nodes = GetAleDiscretization()->GetSolidL2G();
    vector<int> cv;
    for (int i=0;i<s_nodes.size();++i)
      if (i_nodes.find(s_nodes[i])==i_nodes.end())
	{
	  GetGV(gf)(s_nodes[i],0) = 0.0;
	}


    // delete velocity in solid
    int DIM = GetGV(gf).ncomp()==5?2:3;
    assert((DIM==2)||(DIM==3));
    
    for (int i=0;i<s_nodes.size();++i)
      {
    	for (int cc=0;cc<DIM;++cc)
    	  GetGV(gf)(s_nodes[i],cc+1) = 0.0;
      }
    
  }

  /*-------------------------------------------------------*/

  void AleSolver::SetBoundaryVector(VectorInterface& gf) const
  {
    StdSolver::SetBoundaryVector(gf);
    
    const HASHSET<int>&     i_nodes = GetAleDiscretization()->GetInterfaceNodes();
    const vector<int>&      s_nodes = GetAleDiscretization()->GetSolidL2G();
    vector<int> cv;

    for (int i=0;i<s_nodes.size();++i)
      if (i_nodes.find(s_nodes[i])==i_nodes.end())
	{
	  GetGV(gf)(s_nodes[i],0) = 0.0;
	}

    // delete velocity in solid
    int DIM = GetGV(gf).ncomp()==5?2:3;
    for (int i=0;i<s_nodes.size();++i)
      {
    	for (int cc=0;cc<DIM;++cc)
    	  GetGV(gf)(s_nodes[i],cc+1) = 0.0;
      }
    
  }

  void AleSolver::AssembleMatrix(const VectorInterface& gu, double d)
  {
    S2.start();
    StdSolver::AssembleMatrix(gu,d);

    int DIM = GetGV(gu).ncomp()==5?2:3;
    assert((DIM==2)||(DIM==3));
    
    // Modify for pressure zero in Solid-Part
    const HASHSET<int>&     i_nodes = GetAleDiscretization()->GetInterfaceNodes();
    const vector<int>&      s_nodes = GetAleDiscretization()->GetSolidL2G();
    vector<int> cv;
    cv.push_back(0);
    for (int i=0;i<s_nodes.size();++i)
      if (i_nodes.find(s_nodes[i])==i_nodes.end())
	GetMatrix()->dirichlet(s_nodes[i],cv);


    //     In SOlid Part & Interface, 
    //     delete psi-f test-functions & set 
    //     velocity to zero
    cv.clear(); 
    for (int i=0;i<DIM;++i)   cv.push_back(1+i);

    for (HASHSET<int>::const_iterator it = i_nodes.begin();
	 it!=i_nodes.end();++it)
      GetMatrix()->dirichlet(*it,cv);
    for (int i=0;i<s_nodes.size();++i)
      GetMatrix()->dirichlet(s_nodes[i],cv);



    
    
    
    if (__splittingsmoother)
      {
	int ncomp = GetGV(gu).ncomp();
	
	if (ncomp==5)
	  SplitMatrix(dynamic_cast<      SparseBlockMatrix<FMatrixBlock<5> >*> (__AF),
		      dynamic_cast<      SparseBlockMatrix<FMatrixBlock<5> >*> (__AS),
		      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<5> >*> (GetMatrix()));
	else if (ncomp==7)
	  SplitMatrix(//dynamic_cast<SparseBlockMatrix<FluidBlock<3> >*> (__AF),
		      dynamic_cast<SparseBlockMatrix<FMatrixBlock<7> >*> (__AF),
	    
		      dynamic_cast<SparseBlockMatrix<FMatrixBlock<7> >*> (__AS),
		      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<7> >*> (GetMatrix()));
	else abort();

      }
    
	
	
  }


    /*---------------------------------------------------*/

  void AleSolver::AssembleDualMatrix(const VectorInterface& gu, double d)
  {
    __DUALMATRIX=true;
    AssembleMatrix(gu,d);
    __DUALMATRIX=false;


    GetMatrix()->transpose();

    // const HASHSET<int>&     i_nodes = GetAleDiscretization()->GetInterfaceNodes();
    // const vector<int>&      s_nodes = GetAleDiscretization()->GetSolidL2G();
    
    // for (HASHSET<int>::const_iterator it = i_nodes.begin();
    // 	 it!=i_nodes.end();++it)
    //   {
    // 	GetMatrix()->dirichlet_row_col(*it,1,3);
    // 	GetMatrix()->dirichlet_row_col(*it,2,4);
    //   }
    // for (int i=0;i<s_nodes.size();++i)
    //   {
    // 	GetMatrix()->dirichlet_row_col(s_nodes[i],1,3);
    // 	GetMatrix()->dirichlet_row_col(s_nodes[i],2,4);
    //   }
    // DirichletMatrix();
    
    if (__splittingsmoother)
      {
	int ncomp = GetGV(gu).ncomp();
	if (ncomp==5)
	  SplitMatrix(dynamic_cast<      SparseBlockMatrix<FMatrixBlock<5> >*> (__AF),
		      dynamic_cast<      SparseBlockMatrix<FMatrixBlock<5> >*> (__AS),
		      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<5> >*> (GetMatrix()));
	else if (ncomp==7)
	  SplitMatrix(//dynamic_cast<SparseBlockMatrix<FluidBlock<3> >*> (__AF),
		      dynamic_cast<SparseBlockMatrix<FMatrixBlock<7> >*> (__AF),
		      
		      dynamic_cast<SparseBlockMatrix<FMatrixBlock<7> >*> (__AS),
		      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<7> >*> (GetMatrix()));
      }

    ComputeIlu(gu);
  }


  // --------------------------------------------------
  
  void AleSolver::ComputeIlu(const VectorInterface& gu) const
  {
    S5.start();
    if (!__splittingsmoother)
      StdSolver::ComputeIlu(gu);
    S5.stop();
  }
  
  /*-------------------------------------------------------------*/

  MatrixInterface* AleSolver::NewMatrix(int ncomp, const string& matrixtype) 
  {
    if ((matrixtype=="split")&&(!_directsolver))
      {
	__splittingsmoother = true;
	if (ncomp==7)
	  {
	    __AS = new SparseBlockMatrix<FMatrixBlock<7> >;
	    __AF = new SparseBlockMatrix<FMatrixBlock<7> >;
	    //__AF = new SparseBlockMatrix<FluidBlock<3> >;
	    
	    return new SparseBlockMatrix<FMatrixBlock<7> >;
	  }
	else if (ncomp==5)
	  {
	    __AS = new SparseBlockMatrix<FMatrixBlock<5> >;
	    __AF = new SparseBlockMatrix<FMatrixBlock<5> >;
	    return new SparseBlockMatrix<FMatrixBlock<5> >;
	  }
	else abort();
      }
    else return StdSolver::NewMatrix(ncomp,matrixtype);
  }
  
  /*-------------------------------------------------------------*/

  IluInterface* AleSolver::NewIlu(int ncomp, const string& matrixtype) 
  { 
    if ((matrixtype=="split")&&(!_directsolver))
      {
	__splittingsmoother = true;
#ifndef __WITH_UMFPACK__
	cerr << "matrixtype split only with UMFPACK!" << endl;
	abort();
#endif

	if (ncomp==7)
	  {
	    //if (!__splitting_fluid_exact)  __IF = new SparseBlockIlu<FluidBlock<3> >;
	    if (!__splitting_fluid_exact)  __IF = new SparseBlockIlu<FMatrixBlock<7> >;
	    else                           __IF = new SparseUmf<FMatrixBlock<7> > (__AF);

	    if (!__splitting_solid_exact)  __IS = new SparseBlockIlu<FMatrixBlock<7> >;
	    else                           __IS = new SparseUmf<FMatrixBlock<7> > (__AS);
	    
	    return new SparseBlockIlu<FMatrixBlock<7> >;
	  }
	else if (ncomp==5)
	  {
	    if (!__splitting_fluid_exact)  __IF = new SparseBlockIlu<FMatrixBlock<5> >;
	    else                           __IF = new SparseUmf<FMatrixBlock<5> > (__AF);

	    if (!__splitting_solid_exact)  __IS = new SparseBlockIlu<FMatrixBlock<5> >;
	    else                           __IS = new SparseUmf<FMatrixBlock<5> > (__AS);

	    return new SparseBlockIlu<FMatrixBlock<5> >;
	  }
	else assert(0);
      }
    else return StdSolver::NewIlu(ncomp,matrixtype);
  }
  
  
  void AleSolver::Fluidg2l(GlobalVector& Xf, const GlobalVector& x) const
  {
    const std::vector<int>& fluid_l2g = GetAleDiscretization()->GetFluidL2G();
    assert(Xf.ncomp()==x.ncomp());
    assert(Xf.n()==fluid_l2g.size());
    
    for (int i=0;i<fluid_l2g.size();++i) Xf.equ_node(i, fluid_l2g[i], x);
 }
  
  void AleSolver::Solidg2l(GlobalVector& Xs, const GlobalVector& x) const
  {
    const std::vector<int>& solid_l2g = GetAleDiscretization()->GetSolidL2G();
    assert(Xs.ncomp()==x.ncomp());
    assert(Xs.n()==solid_l2g.size());

    for (int i=0;i<solid_l2g.size();++i) Xs.equ_node(i, solid_l2g[i], x);
 }
  
  void AleSolver::Fluidl2g(GlobalVector& x, double s, const GlobalVector& Xf) const
  {
    const std::vector<int>& fluid_l2g = GetAleDiscretization()->GetFluidL2G();
    assert(Xf.ncomp()==x.ncomp());
    assert(Xf.n()==fluid_l2g.size());

    for (int i=0;i<fluid_l2g.size();++i) 
      x.add_node(fluid_l2g[i], s, i, Xf);
 }
  
  void AleSolver::Solidl2g(GlobalVector& x, double s, const GlobalVector& Xs) const
  {
    const std::vector<int>& solid_l2g = GetAleDiscretization()->GetSolidL2G();
    assert(Xs.ncomp()==x.ncomp());
    assert(Xs.n()==solid_l2g.size());
    
    for (int i=0;i<solid_l2g.size();++i) 
      x.add_node(solid_l2g[i], s, i, Xs);
  }
  
  void AleSolver::FluidInterfaceZero(GlobalVector& X, int c) const
  {
    const std::vector<int>& fluid_l2g = GetAleDiscretization()->GetFluidL2G();
    const HASHSET<int>& interface_nodes = GetAleDiscretization()->GetInterfaceNodes();
    
    assert(X.n()==fluid_l2g.size());
    assert(c<X.ncomp());
    
    for (int i=0;i<fluid_l2g.size();++i)
      {
	int n = fluid_l2g[i];
	if (interface_nodes.find(n)!=interface_nodes.end())
	  X(i,c) = 0.0;
      }
  }

  void AleSolver::SolidInterfaceZero(GlobalVector& X, int c) const
  {
    const std::vector<int>& solid_l2g = GetAleDiscretization()->GetSolidL2G();
    const HASHSET<int>& interface_nodes = GetAleDiscretization()->GetInterfaceNodes();
    
    assert(X.n()==solid_l2g.size());
    assert(c<X.ncomp());
    
    for (int i=0;i<solid_l2g.size();++i)
      {
	int n = solid_l2g[i];
	if (interface_nodes.find(n)!=interface_nodes.end())
	  X(i,c) = 0.0;
      }
  }

  /*-------------------------------------------------------*/

  void AleSolver::SolveFluidExact(VectorInterface& x,const VectorInterface& h) const
  {
    assert(__splitting_fluid_exact);
    
    Fluidg2l(__Xf,GetGV(h));         // restrict to fluid-nodes
    
    // for (int c=1;c<__Xf.ncomp();++c) // set v&u rhs to zero on interface
    //   FluidInterfaceZero(__Xf,c);

    ///// set v & u rhs to Dirichlet-Values
    const std::vector<int>& fluid_l2g   = GetAleDiscretization()->GetFluidL2G();
    const HASHSET<int>& interface_nodes = GetAleDiscretization()->GetInterfaceNodes();
    for (int i=0;i<fluid_l2g.size();++i)
      {
	int n = fluid_l2g[i];
	if (interface_nodes.find(n)!=interface_nodes.end())
	  for (int c=1;c<__Xf.ncomp();++c)
	    __Xf(i,c) = GetGV(x)(n,c);
      }

    assert(__IF);                    // solve
    __IF->solve(__Xf);

    for (int c=1;c<__Xf.ncomp();++c) // set v&u rhs to zero on interface
      FluidInterfaceZero(__Xf,c);

    Fluidl2g(GetGV(x),1.0, __Xf);    // update
    
  }  

  /*-------------------------------------------------------*/
  
  void AleSolver::SolveSolidExact(VectorInterface& x,const VectorInterface& h) const
  {
    
    assert(__splitting_solid_exact);
    
    Solidg2l(__Xs,GetGV(h));         // restrict to fluid-nodes
    
    assert(__IS);                    // solve
    __IS->solve(__Xs);

    Solidl2g(GetGV(x),1.0, __Xs);    // update
  }  

  /*-------------------------------------------------------*/

  void AleSolver::SolveFluidIterative(int niter, double TOL, VectorInterface& x,const VectorInterface& h) const
  {
    std::cout << "F\t";
    
    // assert(!__splitting_fluid_exact);
    
    // Fluidg2l(__Xf,GetGV(h));         // restrict to fluid-nodes
    
    // for (int c=1;c<__Xf.ncomp();++c) // set v&u rhs to zero on interface
    //   FluidInterfaceZero(__Xf,c);
    
    // assert(__IF);                    // solve
    // __IF->solve(__Xf);

    // Fluidl2g(GetGV(x),1.0, __Xf);    // update



    // return;
    

    assert(!__splitting_fluid_exact);
    
    Fluidg2l(__Hf,GetGV(h));         // restrict to fluid-nodes
    
    for (int c=1;c<__Hf.ncomp();++c) // set v&u rhs to zero on interface
      FluidInterfaceZero(__Hf,c);
    
    assert(__AF);
    assert(__IF);
    Richardson(*__AF, __Xf, __Hf, *__IF, niter, TOL);
    //BiCGSTAB(*__AF, __Xf, __Hf, *__IF, niter, TOL);
    
    Fluidl2g(GetGV(x),1.0, __Xf);    // update
  }  

  /*-------------------------------------------------------*/
  
  void AleSolver::SolveSolidIterative(int niter, double TOL, VectorInterface& x,const VectorInterface& h) const
  {
    
    assert(!__splitting_solid_exact);
    
    
    Solidg2l(__Hs,GetGV(h));         // restrict to fluid-nodes
    
    assert(__AS);
    assert(__IS);
    Richardson(*__AS, __Xs, __Hs, *__IS, niter, TOL);

    Solidl2g(GetGV(x),1.0, __Xs);    // update
  }  

  /*-------------------------------------------------------*/

  void AleSolver::smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const
  {
    S4.start();
    abort();
    
    double omega = GetSolverData().GetOmega();
    if (__splittingsmoother)
      {

	const std::vector<int>& fluid_l2g = GetAleDiscretization()->GetFluidL2G();
	const std::vector<int>& solid_l2g = GetAleDiscretization()->GetSolidL2G();
	const HASHSET<int>& interface_nodes = GetAleDiscretization()->GetInterfaceNodes();

	const GlobalVector& H = GetGV(h);

	int ncomp = H.ncomp();
	assert((ncomp==5)||(ncomp==7));
	int dim = (ncomp-1)/2;
	assert((dim==2)||(dim==3));

	GlobalVector HF(ncomp, fluid_l2g.size());
	GlobalVector HS(ncomp, solid_l2g.size());

	GlobalVector WF(ncomp, fluid_l2g.size());
	GlobalVector WS(ncomp, solid_l2g.size());

	GlobalVector XF(ncomp, fluid_l2g.size());
	GlobalVector XS(ncomp, solid_l2g.size());

	assert(__IF);
	assert(__IS);

	for(int iter=0; iter<niter; iter++)
	  {
	    
	    ///////////////////////////// FLUID
	    MatrixResidual(h,x,y);
	    
	    for (int i=0;i<fluid_l2g.size();++i) HF.equ_node(i, fluid_l2g[i], H);
	    
	    for (int i=0;i<fluid_l2g.size();++i)
	      {
		int n = fluid_l2g[i];
		if (interface_nodes.find(n)!=interface_nodes.end())
		  for (int ii=0;ii<2*dim;++ii)
		    HF(i,ii+1) = 0.0;
	      }
	    
	    
	    // BICGSTAB

	    double bicgstab_tol = 0.01;
	    int    bicgstab_maxiter = 4;
	    if (__splitting_fluid_exact)
	      {
		XF = HF;
		__IF->solve(XF);
	      }
	    else
	      //	      BiCGSTAB(*__AF, XF, HF, *__IF, bicgstab_maxiter, bicgstab_tol);
	    Richardson(*__AF,XF, HF, *__IF, bicgstab_maxiter, bicgstab_tol);

	    
	    // >>>>>>>>>>>> RICHARDSON!!!
	    /*
	    XF.zero();
	    for (int ii=0;ii<5;++ii)
	      {
		WF = HF;
		if (ii>0)
		  __AF->vmult(WF,XF,-1.0);
		//		cout << "fluid: " << ii << "\t" << WF.norm() << endl;
		__IF->solve(WF);
		XF.add(0.5,WF); 
	      }
	    */
	    // <<<<<<<<<<<< RICHARDSON


	    
	    // for (int i=0;i<fluid_l2g.size();++i) 
	    //   if (interface_nodes.find(fluid_l2g[i])!=interface_nodes.end())
	    //   	HF.scale_node(i,0.5);
	    for (int i=0;i<fluid_l2g.size();++i) 
	      //	      GetGV(x).add_node(fluid_l2g[i], omega, i, HF);
	      GetGV(x).add_node(fluid_l2g[i], omega, i, XF);


	    MatrixResidual(h,x,y);

	    ///////////////////////////// SOLID
	    if (solid_l2g.size()>0)
	      {
		//MatrixResidual(h,x,y);
		for (int i=0;i<solid_l2g.size();++i) HS.equ_node(i, solid_l2g[i], H);
		
		//	    __IS->solve(HS);
		
		bicgstab_tol = 0.001;
		bicgstab_maxiter = 10;
		
		if (__splitting_solid_exact)
		  {
		    XS = HS;
		    __IS->solve(XS);
		  }
		else 
		  //  BiCGSTAB(*__AS, XS, HS, *__IS, bicgstab_maxiter, bicgstab_tol);
		  Richardson(*__AS, XS, HS, *__IS, bicgstab_maxiter, bicgstab_tol);
		
		
		
		// for (int i=0;i<solid_l2g.size();++i) 
		//   if (interface_nodes.find(solid_l2g[i])!=interface_nodes.end())
		// 	HS.scale_node(i,0.5);
		for (int i=0;i<solid_l2g.size();++i) 
		  //GetGV(x).add_node(solid_l2g[i], omega, i, HS);
		  GetGV(x).add_node(solid_l2g[i], omega, i, XS);
	      }
	    
	  }
      }
    else StdSolver::smooth(niter,x,y,h);
    
    S4.stop();
  }
  
  // --------------------------------------------------

  void AleSolver::Visu(const string& name, const VectorInterface& gu, int ii) const
  {
    StdSolver::Visu(name,gu,ii);
    return;
    
    
    // Additional variable which is
    // -1 fluid
    //  0 interface
    //  1 solid

    Chi chi;

    std::string __solid_type;
    DataFormatHandler DFH;
    DFH.insert("solid_type",&__solid_type);
    FileScanner FS(DFH,_paramfile,"Equation");

    chi.BasicInit(__solid_type);
    const GlobalVector& U= GetGV(gu);
    GlobalVector UN;
    UN.ncomp() = U.ncomp()+1;
    UN.resize(U.n());
    for (int i=0;i<U.n();++i)
      {
	for (int c=0;c<U.ncomp();++c)
	  UN(i,c) = U(i,c);

	if (U.ncomp()==2)
	  {
	    const Vertex2d& v=dynamic_cast<const GascoigneMesh2d*> (GetMesh())->vertex2d(i);
	    if (v.x()>=1.0) UN(i,1) = 0.0;
	  }
		
	if (U.ncomp()==5)
	  {
	    const Vertex2d& v=dynamic_cast<const GascoigneMesh2d*> (GetMesh())->vertex2d(i);
	    UN(i,5) = chi(v);
	  }
	else if (U.ncomp()==7)
	  {
	    const Vertex3d& v=dynamic_cast<const GascoigneMesh3d*> (GetMesh())->vertex3d(i);
	    UN(i,7) = static_cast<double> (chi(v));
	  }
      }
    StdSolver::PointVisu(name,UN,ii);
    
    
    

// //     const GlobalVector& U = GetGV(gu);
// //     const GascoigneMesh2d* CM = dynamic_cast<const GascoigneMesh2d*> (GetMesh());
// //     assert(CM);
    
// //     GascoigneMesh2d* M = const_cast<GascoigneMesh2d*> (CM);
// //     assert(M);
// //    std::vector<Vertex2d>& VV = M->GetVertexVector();

// //     for (int i=0;i<M->nnodes();++i)
// //       {
// // 	VV[i][0] += U(i,3);
// // 	VV[i][1] += U(i,4);
// //       }
//      StdSolver::Visu(name,gu,ii);
// //     for (int i=0;i<M->nnodes();++i)
// //       {
// // 	VV[i][0] -= U(i,3);
// // 	VV[i][1] -= U(i,4);
// //       }    
  }
  

  DiscretizationInterface* AleSolver::NewDiscretization(int dimension, const string& discname)
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
	// else if (discname=="AleQ2Lps")            return new AleQ2Lps3d;
	// else
	  return StdSolver::NewDiscretization(dimension, discname);
      }
    else abort();
  }

  void AleSolver::reinit_element_2d(int en, const nvector<int>& indices, 
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

  void AleSolver::reinit_element_3d(int en, const nvector<int>& indices, 
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
    
    
    const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (GetMesh());
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
  

  void AleSolver::ReInitInterface(AleBaseDiscretization* ALEDISC)
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
	    reinit_element_2d(c, M->IndicesOfCell(c), solid_interface_cells, 
			      fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	
	else if ((GetDiscretization()->GetName()=="Q2 Ale 2d Lps")||
		 (GetDiscretization()->GetName()=="Q2 Ale 2d"))
	  for (int c=0;c<M->npatches();++c)
	    reinit_element_2d(c, *(M->IndicesOfPatch(c)), solid_interface_cells, 
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
	    reinit_element_3d(c, M->IndicesOfCell(c), solid_interface_cells, fluid_interface_cells, fluid_cells, solid_cells, 
			      interface_nodes, fluid_nodes, solid_nodes);
	else if (GetDiscretization()->GetName()=="Q2 Ale 3d Lps")
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


    // RESCALE VECTORS
    int ncomp = 1+2*dim;
    __Xf.ncomp () = ncomp;
    __Xf.resize(fluid_l2g.size());
    __Xs.ncomp () = ncomp;
    __Xs.resize(solid_l2g.size());
    __Hf.ncomp () = ncomp;
    __Hf.resize(fluid_l2g.size());
    __Hs.ncomp () = ncomp;
    __Hs.resize(solid_l2g.size());    
  }

  // --------------------------------------------------

  void AleSolver::NewMesh(const MeshInterface* mp)
  {
    cond_computed = false;
    
    //    __splittingsmoother = 0;
    
    // if (mp->nnodes()<1000) __splitting_fluid_exact = true;
    // else __splitting_fluid_exact = false;

    StdSolver::NewMesh(mp);
    ReInitInterface(GetAleDiscretization());
  }
  

#define NCC 5
  void AleSolver::SwapUV(MatrixInterface& M)
  {
    SparseBlockMatrix<FMatrixBlock<NCC> >* SBM =
      dynamic_cast<SparseBlockMatrix<FMatrixBlock<NCC> >* >(&M);
    if (!SBM) return;

    for (int i=0;i<SBM->size();++i)
      {
	FMatrixBlock<NCC> X = *SBM->mat(i);
	for (int j=0;j<2;++j)
	  for (int c=0;c<NCC;++c)
	    {
	      (*SBM->mat(i))(j+1,c) = X(j+3,c);
	      (*SBM->mat(i))(j+3,c) = X(j+1,c);
	    }
      }
    
  }
  

  double AleSolver::LargestEV(MatrixInterface& M)
  {
    const SparseBlockMatrix<FMatrixBlock<NCC> >* SBM =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<NCC> >* >(&M);
    if (!SBM) return 0.0;
    
    const ColumnStencil* ST = 
      dynamic_cast<const ColumnStencil*> (SBM->GetStencil());
    assert(ST);
    
    int n  = ST->n();
    int nc = GetProblemDescriptor()->GetEquation()->GetNcomp();

    if (n<1) return 0.0;    
    double ev = 0.0;

    GlobalVector X,Y;

    X.ncomp()=nc;
    Y.ncomp()=nc;
    X.resize(n);
    Y.resize(n);
    
    X=1.0;
    
    double last_ev;
    
    for (int iter=0;iter<1000;++iter)
      {
	last_ev = ev;
	Y.zero();
	X.zero_comp(0);
	X.zero_comp(1);
	X.zero_comp(2);

	M.vmult(Y,X,1.0);
	Y.zero_comp(0);
	Y.zero_comp(1);
	Y.zero_comp(2);

	ev = Y.norm();
	if ((fabs(last_ev-ev)/ev)<1.e-3) break;
	X=Y;
	X *= 1./ev;
      }
    
    return ev;
  }

  double AleSolver::SmallestEV(MatrixInterface& M)
  {
    const SparseBlockMatrix<FMatrixBlock<NCC> >* SBM =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<NCC> >* >(&M);
    if (!SBM) return 1.0;
    
    const ColumnStencil* ST = 
      dynamic_cast<const ColumnStencil*> (SBM->GetStencil());
    assert(ST);

    int n  = ST->n();
    int nc = GetProblemDescriptor()->GetEquation()->GetNcomp();

    if (n<1) return 1.0;

    double ev = 0.0;

    SparseUmf<FMatrixBlock<NCC> > SU(SBM);
    vector<int> P;
    SU.ConstructStructure(P,*SBM);
    SU.zero();
    SU.copy_entries(SBM);
    SU.compute_ilu();
    

    GlobalVector X,Y;

    X.ncomp()=nc;
    Y.ncomp()=nc;
    X.resize(n);
    Y.resize(n);
    
    X=1.0;

    
    double last_ev;
    
    for (int iter=0;iter<1000;++iter)
      {
	last_ev = ev;
	X.zero_comp(0);
	X.zero_comp(1);
	X.zero_comp(2);
	SU.solve(X);
	X.zero_comp(0);
	X.zero_comp(1);
	X.zero_comp(2);
	ev = X.norm();
	
	if ((fabs(last_ev-ev)/ev)<1.e-3) break;
	X *= 1./ev;
      }
    
    return 1./ev;
  }




  void AleSolver::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
  {
    StdSolver::ConstructInterpolator(I,MT);
    assert(dynamic_cast<AleInterpolator*> (I));
    dynamic_cast<AleInterpolator*> (I)->
      SetInterface(&GetAleDiscretization()->GetInterfaceNodes(),
		   &GetAleDiscretization()->GetFluidG2L(),
		   &GetAleDiscretization()->GetSolidG2L());
  }

  
  


}
