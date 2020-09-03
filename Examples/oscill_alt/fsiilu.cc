#include "fsiilu.h"
#include "fsimatrix.h"
#include "stopwatch.h"
#include "diagblock.h"

namespace Gascoigne {

extern StopWatch _sm_ns, _sm_solid, _sm_ext;

template <int DIM>
FSIIlu<DIM>::FSIIlu(const MatrixInterface* M) {
  __FSIMATRIX= dynamic_cast<const FSIMatrix<DIM>*>(M);
  assert(__FSIMATRIX);

  //    __SM.SetMatrix     (&__FSIMATRIX->GetS11());
  __SES.SetMatrix(&__FSIMATRIX->GetS22());
  //    __F_NS.SetMatrix  (&__FSIMATRIX->GetF_NS());
  //__F_EXT.SetMatrix (&__FSIMATRIX->GetF_EXT());
}

template <int DIM>
void FSIIlu<DIM>::ReInit(const SparseStructureInterface* S) {
}

template <int DIM>
void FSIIlu<DIM>::ConstructStructure(const IntVector&       permfluid,
                                     const IntVector&       permsolid,
                                     const MatrixInterface& __XXX) {
  const FSIMatrix<DIM>* SM= dynamic_cast<const FSIMatrix<DIM>*>(&__XXX);
  assert(SM);

  __SM.ConstructStructure(permsolid, SM->GetS11());

  __SES.ConstructStructure(permsolid, SM->GetS22());

  // // VANKA smoother. We need the patch-structure restricted to the solid:
  // const GascoigneMesh* GM = dynamic_cast<const GascoigneMesh*> (__mesh);
  // assert(GM);
  // nvector<IntVector> __pi;
  // const nvector<IntVector>& pim = GM->GetPatchIndexHandler().GetIndex();
  // for (int p=0;p<pim.size();++p)
  //   {
  // 	IntVector onepatch;
  // 	int max = pim[p].size();
  // 	int j=0;
  // 	for (j=0;j<max;++j)
  // 	  {
  // 	    int sn = Sg2l(pim[p][j]);
  // 	    if (sn==-1) break;
  // 	    onepatch.push_back(sn);
  // 	  }
  // 	if (j==max) __pi.push_back(onepatch);
  //   }
  // __SES.SetMesh(__pi);

  __F_NS.ConstructStructure(permfluid, SM->GetF_NS());
  __F_EXT.ConstructStructure(permfluid, SM->GetF_EXT());

  _xs1.ncomp()= DIM;
  _xs1.resize(SM->GetS11().n());
  _ys1.ncomp()= DIM;
  _ys1.resize(SM->GetS11().n());
  _hs1.ncomp()= DIM;
  _hs1.resize(SM->GetS11().n());

  _xs2.ncomp()= DIM;
  _xs2.resize(SM->GetS22().n());
  _ys2.ncomp()= DIM;
  _ys2.resize(SM->GetS22().n());
  _hs2.ncomp()= DIM;
  _hs2.resize(SM->GetS22().n());

  _xf_ns.ncomp()= DIM + 1;
  _xf_ns.resize(SM->GetF_NS().n());
  _hf_ns.ncomp()= DIM + 1;
  _hf_ns.resize(SM->GetF_NS().n());
  _yf_ns.ncomp()= DIM + 1;
  _yf_ns.resize(SM->GetF_NS().n());
  _xf_ext.ncomp()= DIM;
  _xf_ext.resize(SM->GetF_EXT().n());
  _yf_ext.ncomp()= DIM;
  _yf_ext.resize(SM->GetF_EXT().n());
  _hf_ext.ncomp()= DIM;
  _hf_ext.resize(SM->GetF_EXT().n());
}

template <int DIM>
void FSIIlu<DIM>::copy_entries(const MatrixInterface* A) {
  const FSIMatrix<DIM>* SM= dynamic_cast<const FSIMatrix<DIM>*>(A);
  assert(SM);

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////  SOLID
  //////////////////////////////////////////////////

  const SparseBlockMatrix<DiagBlock<DIM>>& s11= SM->GetS11();
  //    const SparseBlockMatrix<FMatrixBlock<DIM> >&  s12 = SM->GetS12();
  const SparseBlockMatrix<DiagBlock<DIM>>& s21= SM->GetS21();
  //    const SparseBlockMatrix<FMatrixBlock<DIM> >&  s22 = SM->GetS22();

  const ColumnDiagStencil* SS= dynamic_cast<const ColumnDiagStencil*>(s11.GetStencil());
  for (int r= 0; r < SS->n(); ++r) {
    for (int p= SS->start(r); p < SS->stop(r); ++p) {
      const DiagBlock<DIM>& b11= *(s11.mat(p));
      const DiagBlock<DIM>& b21= *(s21.mat(p));

      for (int c1= 0; c1 < DIM; ++c1) {
        double dx= fabs(b11(c1, c1) + b21(c1, c1));
        double sc= b11(c1, c1);
        if (sc != 0)
          assert((fabs(dx / sc) < 1.e-7) || ((b11(c1, c1) == 1) && (b21(c1, c1) == 0)));
      }
    }
  }

  __SM.copy_entries(&SM->GetS11());       // mass matrix
  __SES.copy_entries(&SM->GetS22());      // Elastic 1
  __SES.add_entries(1.0, &SM->GetS12());  // Elastic 2
  //__SES.copy_entries(&SM->GetS12());   // Elastic 2

  // const SparseBlockMatrix<FMatrixBlock<2*DIM > >* SB = &SM->GetS();
  // const ColumnDiagStencil* SS =
  //   dynamic_cast<const ColumnDiagStencil*> (SB->GetStencil());

  // for (int r=0;r<SS->n();++r)
  //   {
  // 	for (int p=SS->start(r);p<SS->stop(r);++p)
  // 	  {
  // 	    int c = SS->col(p);
  // 	    const FMatrixBlock<2*DIM>& B = *(SB->mat(p));
  // 	    for (int c1=0;c1<DIM;++c1)
  // 	      for (int c2=0;c2<DIM;++c2)
  // 		{
  // 		  if (c1!=c2) assert(B(c1,c2)==0);
  // 		  if (c1!=c2) assert(B(c1+DIM,c2)==0);
  // 		  if (c1!=c2) assert(B(c1+DIM,c2+DIM)==0);
  // 		  if (fabs(B(c1,c2)+B(c1+DIM,c2))>1.e-14)
  // 		    assert((B(c1,c2)==1.0)&&(B(c1+DIM,c2)==0.0));
  // 		}
  // 	  }
  //   }

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////  FLUID NS
  //////////////////////////////////////////////////

  __F_NS.copy_entries(&SM->GetF_NS());
  // // modify fluid-NS matrix for Dirichlet on Interface
  // if (1) // THIS IS FOR A DIRECT SOLVER!
  //   {
  // 	nvector<double>& AxF = __F_NS.GetRaw();
  // 	vector<long>&    ApF = __F_NS.GetRawPosition();
  // 	vector<long>&    AcF = __F_NS.GetRawColumn();

  // 	for (HASHSET<int>::const_iterator it = __interface_nodes->begin();
  // 	     it!=__interface_nodes->end();++it)
  // 	  {
  // 	    int global_node = *it;
  // 	    int fluid_node = Fg2l(global_node);
  // 	    assert(fluid_node>=0);

  // 	    // dirichlet for velocity in main equation
  // 	    for (int c=1;c<=DIM;++c)
  // 	      {
  // 		int row = fluid_node*(DIM+1)+c;
  // 		for (int p=ApF[row];p<ApF[row+1];++p)
  // 		  {
  // 		    int col = AcF[p];
  // 		    if (row==col) AxF[p]=1.0;
  // 		    else AxF[p]=0.0;
  // 		  }
  // 	      }
  // 	  }
  //   }

  // modify fluid-NS matrix for Dirichlet on Interface
  if (1) {
    const nvector<int>&      Q= __F_NS.GetQ();
    const ColumnDiagStencil* ST=
      dynamic_cast<const ColumnDiagStencil*>(__F_NS.GetStencil());
    assert(ST);

    for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
         it != __interface_nodes->end();
         ++it) {
      int global_node= *it;
      int fluid_node = Fg2l(global_node);
      assert(fluid_node >= 0);
      assert(fluid_node < Q.size());
      int         row= Q[fluid_node];
      vector<int> cv;
      for (int c= 1; c <= DIM; ++c)
        cv.push_back(c);
      __F_NS.dirichlet(row, cv);
    }
  }

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////  FLUID EXT
  //////////////////////////////////////////////////

  __F_EXT.copy_entries(&SM->GetF_EXT());

  // modify fluid-EXT matrix for Dirichlet on Interface
  // THIS IS FOR SparseBlock
  if (1) {
    const nvector<int>&      Q= __F_EXT.GetQ();
    const ColumnDiagStencil* ST=
      dynamic_cast<const ColumnDiagStencil*>(__F_EXT.GetStencil());
    assert(ST);

    for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
         it != __interface_nodes->end();
         ++it) {
      int global_node= *it;
      int fluid_node = Fg2l(global_node);
      assert(fluid_node >= 0);
      assert(fluid_node < Q.size());
      int         row= Q[fluid_node];
      vector<int> cv;
      for (int c= 0; c < DIM; ++c)
        cv.push_back(c);
      __F_EXT.dirichlet(row, cv);
    }
  }

  // modify fluid-EXT matrix for Dirichlet on Interface
  // if (1)
  //   {
  // 	nvector<double>& AxF = __F_EXT.GetRaw();
  // 	vector<long>&    ApF = __F_EXT.GetRawPosition();
  // 	vector<long>&    AcF = __F_EXT.GetRawColumn();

  // 	for (HASHSET<int>::const_iterator it = __interface_nodes->begin();
  // 	     it!=__interface_nodes->end();++it)
  // 	  {
  // 	    int global_node = *it;
  // 	    int fluid_node = Fg2l(global_node);
  // 	    assert(fluid_node>=0);

  // 	    dirichlet for velocity in main equation
  // 	    for (int c=0;c<DIM;++c)
  // 	      {
  // 		int row = fluid_node*DIM+c;
  // 		for (int p=ApF[row];p<ApF[row+1];++p)
  // 		  {
  // 		    int col = AcF[p];
  // 		    if (row==col) AxF[p]=1.0;
  // 		    else AxF[p]=0.0;
  // 		  }
  // 	      }
  // 	  }
  //   }

  // cout << "Condition Numbers: n = " << __interface_nodes->size() << endl;
  // //    cout << "F:   " << __F.condition_number() << endl;
  // cout << "S:   " << __S.condition_number() << endl;
  // cout << "NS:  " << __F_NS.condition_number() << endl;
  // cout << "EXT: " << __F_EXT.condition_number() << endl << endl;
}

template <int DIM>
void FSIIlu<DIM>::compute_ilu() {
  StopWatch S;
  cout << "ComputIlu: " << __nnodes << endl;
  S.start();
  __F_NS.compute_ilu();
  S.stop();
  cout << "NS: " << S.read() << endl;
  S.reset();
  S.start();
  __F_EXT.compute_ilu();
  S.stop();
  cout << "EX: " << S.read() << endl;
  S.reset();
  S.start();
  __SM.compute_ilu();
  S.stop();
  cout << "SM: " << S.read() << endl;
  S.reset();
  S.start();
  __SES.compute_ilu();
  S.stop();
  cout << "SE: " << S.read() << endl;
}

template <int DIM>
void FSIIlu<DIM>::Solid1_g2l_set(GlobalVector& xs, const GlobalVector& x) const {
  assert(xs.n() == __solid_l2g->size());
  assert(xs.ncomp() == DIM);
  xs.zero();

  int l= 0;
  for (vector<int>::const_iterator it= __solid_l2g->begin(); it != __solid_l2g->end();
       ++it, ++l)
    for (int c= 0; c < DIM; ++c)
      xs(l, c)= x(*it, c + 1);
}
template <int DIM>
void FSIIlu<DIM>::Solid2_g2l_set(GlobalVector& xs, const GlobalVector& x) const {
  assert(xs.n() == __solid_l2g->size());
  assert(xs.ncomp() == DIM);
  xs.zero();

  int l= 0;
  for (vector<int>::const_iterator it= __solid_l2g->begin(); it != __solid_l2g->end();
       ++it, ++l)
    for (int c= 0; c < DIM; ++c)
      xs(l, c)= x(*it, c + 1 + DIM);
}

template <int DIM>
void FSIIlu<DIM>::Solid1_l2g_add(GlobalVector&       x,
                                 const GlobalVector& xs,
                                 double              s) const {
  assert(xs.n() == __solid_l2g->size());
  assert(xs.ncomp() == DIM);
  int l= 0;
  for (vector<int>::const_iterator it= __solid_l2g->begin(); it != __solid_l2g->end();
       ++it, ++l)
    for (int c= 0; c < DIM; ++c)
      x(*it, c + 1)+= s * xs(l, c);
}
template <int DIM>
void FSIIlu<DIM>::Solid2_l2g_add(GlobalVector&       x,
                                 const GlobalVector& xs,
                                 double              s) const {
  assert(xs.n() == __solid_l2g->size());
  assert(xs.ncomp() == DIM);
  int l= 0;
  for (vector<int>::const_iterator it= __solid_l2g->begin(); it != __solid_l2g->end();
       ++it, ++l)
    for (int c= 0; c < DIM; ++c)
      x(*it, c + 1 + DIM)+= s * xs(l, c);
}

template <int DIM>
void FSIIlu<DIM>::Fluid_NS_g2l_set(GlobalVector& xf, const GlobalVector& x) const {
  assert(xf.n() == __fluid_l2g->size());
  assert(xf.ncomp() == DIM + 1);
  int l= 0;
  for (vector<int>::const_iterator it= __fluid_l2g->begin(); it != __fluid_l2g->end();
       ++it, ++l)
    for (int c= 0; c < DIM + 1; ++c)
      xf(l, c)= x(*it, c);
}
template <int DIM>
void FSIIlu<DIM>::Fluid_NS_l2g_add(GlobalVector&       x,
                                   const GlobalVector& xf,
                                   double              s) const {
  assert(xf.n() == __fluid_l2g->size());
  assert(xf.ncomp() == DIM + 1);
  int l= 0;
  for (vector<int>::const_iterator it= __fluid_l2g->begin(); it != __fluid_l2g->end();
       ++it, ++l)
    for (int c= 0; c < DIM + 1; ++c)
      x(*it, c)+= s * xf(l, c);
}

template <int DIM>
void FSIIlu<DIM>::Fluid_EXT_g2l_set(GlobalVector& xf, const GlobalVector& x) const {
  assert(xf.n() == __fluid_l2g->size());
  assert(xf.ncomp() == DIM);
  assert(xf.ncomp() + DIM + 1 == x.ncomp());
  int l= 0;
  for (vector<int>::const_iterator it= __fluid_l2g->begin(); it != __fluid_l2g->end();
       ++it, ++l)
    for (int c= 0; c < DIM; ++c)
      xf(l, c)= x(*it, c + DIM + 1);
}
template <int DIM>
void FSIIlu<DIM>::Fluid_EXT_l2g_add(GlobalVector&       x,
                                    const GlobalVector& xf,
                                    double              s) const {
  assert(xf.n() == __fluid_l2g->size());
  assert(xf.ncomp() == DIM);
  assert(xf.ncomp() + DIM + 1 == x.ncomp());
  int l= 0;
  for (vector<int>::const_iterator it= __fluid_l2g->begin(); it != __fluid_l2g->end();
       ++it, ++l)
    for (int c= 0; c < DIM; ++c)
      x(*it, c + DIM + 1)+= s * xf(l, c);
}

template <int DIM>
void FSIIlu<DIM>::solve(GlobalVector& r) const {
  abort();
}

template <int DIM>
void FSIIlu<DIM>::Precondition(const IluInterface& M, GlobalVector& x) const {
  for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
       it != __interface_nodes->end();
       ++it) {
    int n= Fg2l(*it);
    assert(n >= 0);
    for (int c= 1; c < x.ncomp(); ++c)
      x(n, c)= 0.0;
  }

  M.solve(x);

  for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
       it != __interface_nodes->end();
       ++it) {
    int n= Fg2l(*it);
    assert(n >= 0);
    for (int c= 1; c < x.ncomp(); ++c)
      x(n, c)= 0.0;
  }
}

template <int DIM>
void FSIIlu<DIM>::vmult_NS_dirichlet(const MatrixInterface& A,
                                     GlobalVector&          y,
                                     GlobalVector&          x,
                                     double                 s) const {
  for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
       it != __interface_nodes->end();
       ++it) {
    int n= Fg2l(*it);
    assert(n >= 0);
    for (int c= 1; c < y.ncomp(); ++c)
      x(n, c)= 0.0;
  }
  A.vmult(y, x, s);
  for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
       it != __interface_nodes->end();
       ++it) {
    int n= Fg2l(*it);
    assert(n >= 0);
    for (int c= 1; c < y.ncomp(); ++c)
      y(n, c)= 0.0;
  }
}

template <int DIM>
void FSIIlu<DIM>::vmult_EXT_dirichlet(const MatrixInterface& A,
                                      GlobalVector&          y,
                                      GlobalVector&          x,
                                      double                 s) const {
  for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
       it != __interface_nodes->end();
       ++it) {
    int n= Fg2l(*it);
    assert(n >= 0);
    x.zero_node(n);
  }
  A.vmult(y, x, s);
  for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
       it != __interface_nodes->end();
       ++it) {
    int n= Fg2l(*it);
    assert(n >= 0);
    y.zero_node(n);
  }
}

template <int DIM>
void FSIIlu<DIM>::Jacobi(const MatrixInterface& A, GlobalVector& x) const {
  return;
  const SparseBlockMatrix<DiagBlock<DIM>>* A1=
    dynamic_cast<const SparseBlockMatrix<DiagBlock<DIM>>*>(&A);
  if (A1) {
    assert(x.ncomp() == DIM);
    const ColumnDiagStencil* ST= dynamic_cast<const ColumnDiagStencil*>(A1->GetStencil());

    assert(x.n() == ST->n());

    for (int r= 0; r < ST->n(); ++r) {
      int                   p= ST->diag(r);
      const DiagBlock<DIM>& B= *(A1->mat(p));
      for (int c= 0; c < x.ncomp(); ++c) {
        if (B(c, c) == 0)
          assert(x(r, c) == 0);
        else
          x(r, c)/= B(c, c);
      }
    }

  } else
    abort();
}

template <int DIM>
int FSIIlu<DIM>::CG(const MatrixInterface& A,
                    GlobalVector&          x,
                    const GlobalVector&    b,
                    const IluInterface&    precond,
                    int&                   max_iter,
                    double&                tol) const {
  x.zero();
  GlobalVector r= b;
  // A.vmult(r,x,-1.0);
  double res = r.norm();
  double goal= res * tol;
  //    cout << "CG: 0\t" << res << endl;

  GlobalVector h= r;
  //    Jacobi(A,h);

  precond.solve(h);

  GlobalVector d= h;
  GlobalVector z= x;

  int iter= 0;
  for (iter= 1; iter <= max_iter; ++iter) {
    z.zero();
    vmult_EXT_dirichlet(A, z, d, 1.0);
    //	A.vmult(z,d,1.0);

    double a1= r * h;
    double a2= d * z;
    if (a1 != 0)
      assert(a2 != 0);

    double alpha= a1 / a2;
    if (a1 == 0)
      alpha= 0.0;
    if (alpha == 0)
      return 1;

    x.add(alpha, d);
    r.add(-alpha, z);
    res= r.norm();
    //	cout << "CG: " << iter << "\t"  << res << endl;
    if (res < goal)
      return 0;
    h= r;
    //	Jacobi(A,h);
    precond.solve(h);

    double b1= r * h;
    double b2= a1;
    assert(b2 != 0);
    double beta= b1 / b2;
    d.sequ(beta, 1, h);
  }
  return 1;
}

template <int DIM>
int FSIIlu<DIM>::CGSolid(const MatrixInterface& A,
                         GlobalVector&          x,
                         const GlobalVector&    b,
                         const IluInterface&    precond,
                         int&                   max_iter,
                         double&                tol) const {
  x.zero();
  GlobalVector r= b;
  // A.vmult(r,x,-1.0);
  double res = r.norm();
  double goal= res * tol;
  //    cout << "CG: 0\t" << res << endl;

  GlobalVector h= r;
  //    Jacobi(A,h);

  precond.solve(h);

  GlobalVector d= h;
  GlobalVector z= x;

  int iter= 0;
  for (iter= 1; iter <= max_iter; ++iter) {
    z.zero();
    A.vmult(z, d, 1.0);

    double a1= r * h;
    double a2= d * z;
    if (a1 != 0)
      assert(a2 != 0);

    double alpha= a1 / a2;
    if (a1 == 0)
      alpha= 0.0;
    if (alpha == 0)
      return 1;

    x.add(alpha, d);
    r.add(-alpha, z);
    res= r.norm();
    if (res < goal)
      return 0;
    h= r;
    //	Jacobi(A,h);
    precond.solve(h);

    double b1= r * h;
    double b2= a1;
    assert(b2 != 0);
    double beta= b1 / b2;
    d.sequ(beta, 1, h);
  }
  return 1;
}

// template<int DIM>
// int  FSIIlu<DIM>::BiGSTABL(int L,
// 			      const MatrixInterface &A,
// 			      GlobalVector &x, const GlobalVector &b,
// 			      const IluInterface &M,
// 			      int &max_iter, double &tol) const
// {
//   int n     = x.n();
//   int ncomp = x.ncomp();

//   vector<GlobalVector> r(L+1,GlobalVector(ncomp,n) );
//   vector<GlobalVector> u(L+1,GlobalVector(ncomp,n) );
//   assert(r[0].n()==n);
//   assert(r[0].ncomp()==ncomp);

//   double bnrm = b.norm();
//   if (bnrm == 0.0) bnrm = 1.0;

//   int iter = 0;

//   vector<double> gamma(L+1), gamma_p(L+1), gamma_pp(L+1), tau(L*L), sigma(L+1);

//   int ierr = 0; // error code to return, if any
//   const double breaktol = 1e-30;

//   // rtilde = r[0] = b - Ax
//   r[0] = b;
//   vmult_NS_dirichlet(A,r[0],x,-1);
//   GlobalVector rtilde = r[0];

//   { /* Sleipjen normalizes rtilde in his code; it seems to help slightly */
//     double s = rtilde.norm();
//     assert(s>0);
//     rtilde *= 1.0/s;
//   }

//   u[0].zero();

//   double rho = 1.0, alpha = 0, omega = 1;

//   double resid;
//   while ((resid = r[0].norm()) > tol * bnrm)
//     {
// 	++iter;
// 	//	cout << "residual " << iter << " " << resid/bnrm << endl;

// 	rho = -omega * rho;
// 	for (int j = 0; j < L; ++j)
// 	  {
// 	    if (fabs(rho) < breaktol)
// 	      {
// 		ierr = -1; goto finish;
// 	      }
// 	    double rho1 = r[j]*rtilde;
// 	    double beta = alpha * rho1 / rho;
// 	    rho = rho1;
// 	    for (int i = 0; i <= j; ++i)
// 	    {
// 	      //for (int m = 0; m < n; ++m) u[i][m] = r[i][m] - beta * u[i][m];
// 	      u[i].sadd(-beta,1.0,r[i]);
// 	    }
// 	    //   A(u[j], u[j+1], Adata);
// 	    u[j+1].zero();
// 	  vmult_NS_dirichlet(A,u[j+1],u[j],1.0);

// 	  alpha = rho / (u[j+1]*rtilde);
// 	  for (int i = 0; i <= j; ++i)
// 	    //xpay(n, r[i], -alpha, u[i+1]);
// 	    r[i].add(-alpha,u[i+1]);

// 	  //A(r[j], r[j+1], Adata);
// 	  r[j+1].zero();
// 	  vmult_NS_dirichlet(A,r[j+1],r[j],1.0);
// 	  //xpay(n, x, alpha, u[0]);
// 	  x.add(alpha,u[0]);
// 	  }

// 	for (int j = 1; j <= L; ++j)
// 	  {
// 	    for (int i = 1; i < j; ++i)
// 	      {
// 		int ij = (j-1)*L + (i-1);
// 		tau[ij] = r[j]*r[i]/sigma[i];
// 		r[j].add(-tau[ij],r[i]);
// 	      }
// 	    sigma[j] = r[j]*r[j];
// 	    gamma_p[j] = r[0]*r[j]/sigma[j];
// 	  }

// 	omega = gamma[L] = gamma_p[L];
// 	for (int j = L-1; j >= 1; --j)
// 	  {
// 	    gamma[j] = gamma_p[j];
// 	    for (int i = j+1; i <= L; ++i)
// 	      gamma[j] -= tau[(i-1)*L + (j-1)] * gamma[i];
// 	  }
// 	for (int j = 1; j < L; ++j)
// 	  {
// 	    gamma_pp[j] = gamma[j+1];
// 	    for (int i = j+1; i < L; ++i)
// 	      gamma_pp[j] += tau[(i-1)*L + (j-1)] * gamma[i+1];
// 	  }
// 	//xpay(n, x, gamma[1], r[0]);
// 	x.add(gamma[1],r[0]);
// 	//xpay(n, r[0], -gamma_p[L], r[L]);
// 	r[0].add(-gamma_p[L],r[L]);
// 	//xpay(n, u[0], -gamma[L], u[L]);
// 	u[0].add(-gamma[L],u[L]);

// 	for (int j = 1; j < L; ++j)
// 	  { /* TODO: use blas DGEMV (for L > 2) */
// 	    //xpay(n, x, gamma_pp[j], r[j]);
// 	    x.add(gamma_pp[j],r[j]);
// 	    //xpay(n, r[0], -gamma_p[j], r[j]);
// 	    r[0].add(-gamma_p[j],r[j]);
// 	    //xpay(n, u[0], -gamma[j], u[j]);
// 	    u[0].add(-gamma[j],u[j]);
// 	  }

// 	if (iter == max_iter) { ierr = 1; break; }
//     }

//   //    cout << "final residual = " << r[0].norm()/bnrm << endl;

// finish:
//   max_iter = iter;
//   return ierr;
// }

template <int DIM>
void FSIIlu<DIM>::solve_fluid(GlobalVector& r) const {
  _sm_ext.start();

  // set zero dirichlet on interface
  for (HASHSET<int>::const_iterator it= __interface_nodes->begin();
       it != __interface_nodes->end();
       ++it)
    for (int c= 1; c < r.ncomp(); ++c)
      r(*it, c)= 0.0;

  // solve extension

  // restrict rhs

  string smooth_ext= "cg";
  if (smooth_ext == "onestep") {
    Fluid_EXT_g2l_set(_xf_ext, r);
    __F_EXT.solve(_xf_ext);
  } else if (smooth_ext == "richardson") {
    Fluid_EXT_g2l_set(_yf_ext, r);
    _xf_ext.zero();
    GlobalVector hf_ext= _xf_ext;
    for (int it= 0; it < 4; ++it) {
      hf_ext= _yf_ext;
      if (it > 0)
        vmult_NS_dirichlet(__FSIMATRIX->GetF_EXT(), hf_ext, _xf_ext, -1.0);
      __F_EXT.solve(hf_ext);
      _xf_ext.add(1.0, hf_ext);
    }
  } else if (smooth_ext == "cg") {
    Fluid_EXT_g2l_set(_yf_ext, r);
    _xf_ext.zero();
    int    maxiter= 12;
    double tol    = 1.e-3;
    maxiter       = CG(__FSIMATRIX->GetF_EXT(), _xf_ext, _yf_ext, __F_EXT, maxiter, tol);
  }

  _sm_ext.stop();

  _sm_ns.start();

  // solve NS
  // NS u_ns = r_ns - ALE _xf_ext = _yf_ns

  string smooth= "richardson";
  if (smooth == "onestep")  // one-step ILU
  {
    Fluid_NS_g2l_set(_xf_ns, r);
    __FSIMATRIX->GetF_ALE().vmult(_xf_ns, _xf_ext, -1.0);
    __F_NS.solve(_xf_ns);
  } else if (smooth == "richardson")  // simple richardson
  {
    Fluid_NS_g2l_set(_yf_ns, r);
    __FSIMATRIX->GetF_ALE().vmult(_yf_ns, _xf_ext, -1.0);
    int maxiter= 20;
    _xf_ns.zero();
    for (int iter= 0; iter < maxiter; ++iter) {
      _hf_ns= _yf_ns;
      vmult_NS_dirichlet(__FSIMATRIX->GetF_NS(), _hf_ns, _xf_ns, -1.0);
      //	    __FSIMATRIX->GetF_NS().vmult(_hf_ns,_xf_ns,-1.0);
      //	    cout << iter << "\t" << _hf_ns.norm() << endl;

      // if (iter>0)
      //   __FSIMATRIX->GetF_NS().vmult(_hf_ns,_xf_ns,-1.0);
      __F_NS.solve(_hf_ns);
      _xf_ns.add(1., _hf_ns);
    }
  } else if (smooth == "bicgstab")  // bicg richardson
  {
    int    maxiter= 10;
    double tol    = 1.e-7;

    _xf_ns.zero();

    Fluid_NS_g2l_set(_yf_ns, r);
    __FSIMATRIX->GetF_ALE().vmult(_yf_ns, _xf_ext, -1.0);
    BiCGSTAB(__FSIMATRIX->GetF_NS(), _xf_ns, _yf_ns, __F_NS, maxiter, tol);
    // int res = BiCGSTABL(8,__FSIMATRIX->GetF_NS(),
    // 		   _xf_ns, _yf_ns, __F_NS,maxiter,tol);

  } else
    abort();
  _sm_ns.stop();

  r.zero();
  Fluid_NS_l2g_add(r, _xf_ns, 1.0);
  Fluid_EXT_l2g_add(r, _xf_ext, 1.0);
}

template <int DIM>
int FSIIlu<DIM>::BiCGSTAB(const MatrixInterface& A,
                          GlobalVector&          x,
                          const GlobalVector&    b,
                          const IluInterface&    M,
                          int&                   max_iter,
                          double&                tol) const {
  int          n    = x.n();
  int          ncomp= x.ncomp();
  double       resid;
  double       rho_1, rho_2, alpha, beta, omega;
  GlobalVector p(ncomp, n), phat(ncomp, n), s(ncomp, n), shat(ncomp, n), t(ncomp, n),
    v(ncomp, n);

  double normb= b.norm();
  if (normb == 0.0)
    normb= 1;

  GlobalVector r= b;

  //    vmult_NS_dirichlet(A,r,x,-1.0);
  A.vmult(r, x, -1.0);

  GlobalVector rtilde= r;

  if ((resid= r.norm() / normb) <= tol) {
    tol     = resid;
    max_iter= 0;
    return 0;
  }
  //    cout << 0 << "\t" << resid << endl;
  for (int i= 1; i <= max_iter; i++) {
    rho_1= rtilde * r;
    if (rho_1 == 0) {
      tol= r.norm() / normb;
      return 2;
    }
    if (i == 1)
      p= r;
    else {
      beta= (rho_1 / rho_2) * (alpha / omega);
      //	  p = r + beta(0) * (p - omega(0) * v);
      p.sadd(beta, 1.0, r);
      p.add(-beta * omega, v);
    }
    phat= p;
    M.solve(phat);
    v.zero();
    //      vmult_NS_dirichlet(A,v,phat,1.0);
    A.vmult(v, phat, 1.0);

    alpha= rho_1 / (rtilde * v);
    s    = r;
    s.add(-alpha, v);

    if ((resid= s.norm() / normb) < tol) {
      x.add(alpha, phat);
      tol= resid;
      return 0;
    }
    shat= s;
    M.solve(shat);
    t.zero();
    //      vmult_NS_dirichlet(A,t,shat,1.0);
    A.vmult(t, shat, 1.0);

    omega= (t * s) / (t * t);
    x.add(alpha, phat, omega, shat);
    r= s;
    r.add(-omega, t);

    rho_2= rho_1;
    if ((resid= r.norm() / normb) < tol) {
      tol     = resid;
      max_iter= i;
      return 0;
    }
    //      cout << i << "\t" << resid << endl;
    if (omega == 0) {
      tol= r.norm() / normb;
      return 3;
    }
  }

  tol= resid;
  return 1;
}

template <int DIM>
int FSIIlu<DIM>::BiCGSTAB(const MatrixInterface& A1,
                          const MatrixInterface& A2,
                          GlobalVector&          x,
                          const GlobalVector&    b,
                          const IluInterface&    M,
                          int&                   max_iter,
                          double&                tol) const {
  int          n    = x.n();
  int          ncomp= x.ncomp();
  double       resid;
  double       rho_1, rho_2, alpha, beta, omega;
  GlobalVector p(ncomp, n), phat(ncomp, n), s(ncomp, n), shat(ncomp, n), t(ncomp, n),
    v(ncomp, n);

  double normb= b.norm();
  if (normb == 0.0)
    normb= 1;

  GlobalVector r= b;
  // for (int i=0;i<__boundarynodes.size();++i)
  //   r.zero_node(__boundarynodes[i]);

  // for (int i=0;i<__boundarynodes.size();++i)
  //   x.zero_node(__boundarynodes[i]);

  //    vmult_NS_dirichlet(A,r,x,-1.0);
  A1.vmult(r, x, -1.0);
  A2.vmult(r, x, -1.0);

  GlobalVector rtilde= r;

  if ((resid= r.norm() / normb) <= tol) {
    tol     = resid;
    max_iter= 0;
    return 0;
  }
  for (int i= 1; i <= max_iter; i++) {
    rho_1= rtilde * r;
    if (rho_1 == 0) {
      tol     = r.norm() / normb;
      max_iter= i;
      return 2;
    }
    if (i == 1)
      p= r;
    else {
      beta= (rho_1 / rho_2) * (alpha / omega);
      //	  p = r + beta(0) * (p - omega(0) * v);
      p.sadd(beta, 1.0, r);
      p.add(-beta * omega, v);
    }
    phat= p;
    M.solve(phat);
    v.zero();
    //      vmult_NS_dirichlet(A,v,phat,1.0);
    // for (int i=0;i<__boundarynodes.size();++i)
    // 	phat.zero_node(__boundarynodes[i]);

    A1.vmult(v, phat, 1.0);
    A2.vmult(v, phat, 1.0);

    alpha= rho_1 / (rtilde * v);
    s    = r;
    s.add(-alpha, v);

    if ((resid= s.norm() / normb) < tol) {
      x.add(alpha, phat);
      tol     = resid;
      max_iter= i;
      return 0;
    }
    shat= s;
    M.solve(shat);
    t.zero();
    //      vmult_NS_dirichlet(A,t,shat,1.0);
    // for (int i=0;i<__boundarynodes.size();++i)
    // 	shat.zero_node(__boundarynodes[i]);

    A1.vmult(t, shat, 1.0);
    A2.vmult(t, shat, 1.0);

    omega= (t * s) / (t * t);
    x.add(alpha, phat, omega, shat);
    r= s;
    r.add(-omega, t);

    rho_2= rho_1;
    if ((resid= r.norm() / normb) < tol) {
      tol     = resid;
      max_iter= i;
      return 0;
    }
    if (omega == 0) {
      tol     = r.norm() / normb;
      max_iter= i;
      return 3;
    }
  }

  tol= resid;

  return 1;
}

template <int DIM>
void FSIIlu<DIM>::mg_smooth_ses(GlobalVector& x,
                                GlobalVector& y,
                                GlobalVector& h,
                                int           iter) const {
  for (int i= 0; i < iter; ++i) {
    h= y;
    __FSIMATRIX->GetS12().vmult(h, x, -1.0);
    __FSIMATRIX->GetS22().vmult(h, x, -1.0);
    __SES.solve(h);
    x.add(1.0, h);
  }
}

template <int DIM>
void FSIIlu<DIM>::mg_restrict_ses(GlobalVector& Y, const GlobalVector& h) const {
  Y.ncomp()= h.ncomp();
  Y.resize(__interpolator.GetC2F().size());
  __interpolator.restrict_zero(Y, h);
}

template <int DIM>
void FSIIlu<DIM>::mg_prolongate_add_ses(GlobalVector& x, const GlobalVector& X) const {
  __interpolator.prolongate_add(x, X);
}

template <int DIM>
void FSIIlu<DIM>::mg_ses(GlobalVector& x, GlobalVector& y, GlobalVector& h) const {
  //    cout << "mg_ses" << thislevel << " " << coarselevel << " " << x.n() << endl;
  // for (int i=0;i<__boundarynodes.size();++i)
  //   y.zero_node(__boundarynodes[i]);

  if (thislevel == coarselevel) {
    cout << y.norm() << endl;
    int    maxiter= 20;
    double tol    = 1.e-4;
    x.zero();
    BiCGSTAB(__FSIMATRIX->GetS22(), __FSIMATRIX->GetS12(), x, y, __SES, maxiter, tol);
  } else {
    assert(thislevel > 0);
    int    maxiter= 4;
    double tol    = 1.e-4;
    x.zero();
    BiCGSTAB(__FSIMATRIX->GetS22(), __FSIMATRIX->GetS12(), x, y, __SES, maxiter, tol);

    // smooth pre
    //	mg_smooth_ses(x,y,h,100); // 4 steps

    // // residual
    h= y;
    __FSIMATRIX->GetS12().vmult(h, x, -1.0);
    __FSIMATRIX->GetS22().vmult(h, x, -1.0);

    // // restrict
    GlobalVector Y, X, H;
    mg_restrict_ses(Y, h);

    X.ncomp()= Y.ncomp();
    X.resize(Y.n());
    X.zero();
    H.ncomp()= Y.ncomp();
    H.resize(Y.n());
    H.zero();

    // // coarse mesh
    assert(__FSI_H[thislevel - 1] != NULL);
    __FSI_H[thislevel - 1]->mg_ses(X, Y, H);

    // // prolongate
    mg_prolongate_add_ses(x, X);
    // for (int i=0;i<__boundarynodes.size();++i)
    //   x.zero_node(__boundarynodes[i]);

    // // smooth post
    // mg_smooth_ses(x,y,h,2);
  }
}

template <int DIM>
void FSIIlu<DIM>::solve_solid(GlobalVector& r) const {
  _sm_solid.start();

  Solid1_g2l_set(_ys1, r);
  Solid2_g2l_set(_ys2, r);

  // first Problem:   (ES + B) xs_2 = ys1+ys2

  string smooth_first= "onestep";
  if (smooth_first == "onestep")  // one-step ILU
  {
    _xs2= _ys2;
    _xs2.add(1.0, _ys1);
    __SES.solve(_xs2);
  } else if (smooth_first == "mg") {
    _ys2.add(1.0, _ys1);
    _xs2.zero();
    mg_ses(_xs2, _ys2, _hs2);
  } else if (smooth_first == "richardson") {
    _ys2.add(1.0, _ys1);
    int maxiter= 50;
    _xs2.zero();
    double tol= 1.e-5;

    double res0= _ys2.norm();
    for (int iter= 0; iter < maxiter; ++iter) {
      _hs2= _ys2;
      if (iter > 0) {
        __FSIMATRIX->GetS22().vmult(_hs2, _xs2, -1.0);
        __FSIMATRIX->GetS12().vmult(_hs2, _xs2, -1.0);
      }
      double res= _hs2.norm();
      if (res <= res0 * tol)
        break;
      __SES.solve(_hs2);
      _xs2.add(1.0, _hs2);
    }
  } else if (smooth_first == "bicgstab") {
    _ys2.add(1.0, _ys1);
    int maxiter= 200;
    _xs2.zero();
    double tol= 1.e-6;

    BiCGSTAB(
      __FSIMATRIX->GetS22(), __FSIMATRIX->GetS12(), _xs2, _ys2, __SES, maxiter, tol);
    //	cout << "XX: " << _xs2.n() << " " << maxiter << " " << tol << endl;
  } else
    abort();

  // second problem M xs_1 = ys1 - ES xs_2
  __FSIMATRIX->GetS12().vmult(_ys1, _xs2, -1.0);

  string smooth_second= "cg";
  if (smooth_second == "cg") {
    int    max= 20;
    double tol= 1.e-5;
    CGSolid(__FSIMATRIX->GetS11(), _xs1, _ys1, __SM, max, tol);
  }

  r.zero();
  Solid1_l2g_add(r, _xs1, 1.0);
  Solid2_l2g_add(r, _xs2, 1.0);

  _sm_solid.stop();
}

template class FSIIlu<2>;
template class FSIIlu<3>;

}  // namespace Gascoigne
