#include  "ale-dual.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
  template<int DIM>
  AleDual<DIM>::AleDual()
  {
    abort();
  }
  
  /*-----------------------------------------*/
  
  template<int DIM>
  AleDual<DIM>::AleDual(const ParamFile* pf) : Equation()
  {
    std::string __solid_type;
    
    __U = 0;
    DataFormatHandler DFH;
    DFH.insert("nu_f", &__nu_f,1.);
    DFH.insert("mu_s", &__mu_s,1.);
    DFH.insert("lambda_s", &__lambda_s,1.);
    DFH.insert("rho_f",&__rho_f,1.);
    DFH.insert("rho_s",&__rho_s,1.);

    DFH.insert("alpha_p",&__alpha_p0,1.);
    DFH.insert("alpha_u",&__alpha_u0,1.);
    DFH.insert("alpha_v",&__alpha_v0,1.);

    DFH.insert("delta_lap_p",&__delta_lap_p0,0.);
    DFH.insert("delta_lps_p",&__delta_lps_p0,0.);
    DFH.insert("delta_lps_v",&__delta_lps_v0,0.);
    DFH.insert("solid_type" ,&__solid_type);
    FileScanner FS(DFH,pf,"Equation");

    __chi.BasicInit(__solid_type);


    __F.resize(DIM,DIM);
    __Ftilde.resize(DIM,DIM);
    __nablaV_Ftilde.resize(DIM,DIM);

    std::cout << "Parameters: " << std::endl;
    std::cout << __rho_f << "\t" << __rho_s << "\t" << __nu_f << "\t" << __mu_s << "\t" << __lambda_s << std::endl;
    std::cout << __alpha_p0 << "\t" << __alpha_u0 << "\t" << __alpha_v0 << std::endl;
    
    
  }
  
  /*-----------------------------------------*/

  template<int DIM>
  void AleDual<DIM>::compute_transformation(const FemFunction& U) const
  {
    assert(DIM==2);
    
    // F
    __F.zero();
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	__F(i,j) = U[i+1+DIM][j+1];     // ok, ok wicki
    for (int i=0;i<DIM;++i) 
      __F(i,i) += 1.0;                  // ok, ok wicki

    // J
    __J = 0.0;
    __J = __F.det();

    // Ftilde
    __Ftilde(0,0) =  __F(1,1);
    __Ftilde(1,1) =  __F(0,0);
    __Ftilde(1,0) = -__F(1,0);
    __Ftilde(0,1) = -__F(0,1);          // ok , ok wicki


    // convektion V in Halb-ALE Coord.  (Ftilde v) \nabla v
    __conv_v=0;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	__conv_v[i] += U[j+1].m() * U[i+1][j+1]; // ok , ok wicki

    // sigma_f : \nabla v * F_tilde --> ist nur der `halbe' Tensor
    __nablaV_Ftilde.zero();
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	for (int k=0;k<DIM;++k)
	  __nablaV_Ftilde(i,j) += U[i+1][k+1] * __Ftilde(k,j); // ok , ok wicki

    // divergence v in ALE: div(JF^{-1} v)
    __divergence_v = 0;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
 	__divergence_v += __Ftilde(i,j) * U[j+1][i+1]; // ok , ok wicki

    // trace E
    __trace_E = -0.5 * DIM;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	__trace_E += 0.5 * __F(j,i)*__F(j,i);  // ok wicki
  }

  template<int DIM>
  double AleDual<DIM>::DU_trace_E(int k, const TestFunction& M) const
  {
    double res = 0.0;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	res += DU_F(i,j,k,M)*__F(i,j);
    return res;
  }

  template<int DIM>
  double AleDual<DIM>::DU_F(int i,int j, int k, const TestFunction& M) const
  {
    if (i!=k) return 0.0;
    return M[j+1];
  }
  

  template<int DIM>
  double AleDual<DIM>::DU_divergence_v(int k, const FemFunction& U, const TestFunction& M) const
  {
    double res = 0.0;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	res += DU_Ftilde(i,j,k,M) * U[j+1][i+1];
    return res;
  }
  template<int DIM>
  double AleDual<DIM>::DV_divergence_v(int k, const TestFunction& M) const
  {
    double res = 0.0;
    for (int i=0;i<DIM;++i)
      res +=  __Ftilde(i,k) * M[i+1];
    return res;
  }
  

  template<int DIM>
  double AleDual<DIM>::DU_J(int k, const FemFunction& U, const TestFunction& M) const
  {
    if (k==0) return M.x() * (1.0+U[1+DIM+1].y()) - M.y() * U[1+DIM+1].x();
    if (k==1) return M.y() * (1.0+U[1+DIM+0].x()) - M.x() * U[1+DIM+0].y();
    abort();
  }
  


  template<int DIM>
  double AleDual<DIM>::DU_nablaV_Ftilde(int i, int j, int k, const FemFunction& U, const TestFunction& M) const
  {
    double res = 0.0;
    for (int l=0;l<DIM;++l)
      res += U[i+1][l+1] * DU_Ftilde(l,j,k,M);
    return res;
  }
  template<int DIM>
  double AleDual<DIM>::DV_nablaV_Ftilde(int i, int j, int k, const TestFunction& M) const 
  {
    if (i!=k) return 0.0;
    double res = 0.0;
    for (int l=0;l<DIM;++l)  res += M[l+1] * __Ftilde(l,j);
    return res;
  }

  // --------------------------------------------------

  template<int DIM>
  double AleDual<DIM>::DV_conv_v(int i, int k, const FemFunction& U, const TestFunction& M) const
  {
    double res = M.m() * U[i+1][k+1];
    if (k==i)
      for (int j=0;j<DIM;++j) 
	res += U[j+1].m() * M[j+1];
    return res;
  }

  template<int DIM>
  double AleDual<DIM>::DU_Ftilde(int i, int j, int k, const TestFunction& M) const
  {
    if ( (i==0)&&(j==0)&&(k==1) ) return  M.y();
    if ( (i==0)&&(j==1)&&(k==0) ) return -M.y();
    if ( (i==1)&&(j==0)&&(k==1) ) return -M.x();
    if ( (i==1)&&(j==1)&&(k==0) ) return  M.x();
    return 0.0;
  }

  /*-----------------------------------------*/

  template<int DIM>
  void AleDual<DIM>::point(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    __domain = __chi(v);


    __alpha_v = __alpha_v0;


//     double dist = 0;
//     if (v.x()>0.6)
//       dist = sqrt( (v.x()-0.6)*(v.x()-0.6) + (v.y()-0.2)*(v.y()-0.2));
//     else if (v.x()<0.25)
//       dist = sqrt( (v.x()-0.25)*(v.x()-0.25) + (v.y()-0.2)*(v.y()-0.2));
//     else
//       dist = fabs(v.y()-0.2);

    double dist = 0.0;
    double dx=0;
    double dy=0;
    if (v.x()<0) dx = fabs(v.x());
    if (v.y()>1.0) dx = fabs(v.x()-1.0);
    if (v.y()>1.0) dy = fabs(v.y()-1.0);
    dist = sqrt(dx*dx+dy*dy);
    
    // 1. Möglichkeit
    //__alpha_u =  __alpha_u0 * h*h;// __alpha_u0*(1.0 + 100.0 * exp(-10.0*dist))  * h*h;
  
    // 2. Möglichkeit
    //__alpha_u =  __alpha_u0*(1.0 + 100.0 * exp(-1.0*dist))  * h*h;

    // 3. Möglichkeit
    __alpha_u =  __alpha_u0;
    
    // 1. Möglichkeit
    //__alpha_p = __alpha_p0 * h*h;//;///h;//;// * h*h;
    
    // 2. Möglichkeit
    //__alpha_p =  __alpha_p0*(1.0 + 100.0 * exp(-1.0*dist))  * h;

    // 3. Möglichkeit
    __alpha_p =  __alpha_p0;



    double vm = 0;
    for (int d=0;d<DIM;++d) vm += U[d+1].m() * U[d+1].m();
    vm = sqrt(vm);
    vm = 0.0;
    
    __delta_lap_p =  __delta_lap_p0/__rho_f / (__nu_f/ (h*h));
    if (__domain>=0) __delta_lap_p = 0;

    compute_transformation(U);
  }

  /*-----------------------------------------*/

  template<int DIM>
  void AleDual<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    assert(__U);
    if (__domain<0) // fluid
      {
	// stab
	for (int j=0;j<DIM;++j)
	  b[0] += __delta_lap_p * U[0][j+1] * N[j+1];

	// new:	
	// dual stress term including pressure
	for (int i=0;i<DIM;++i)	  
	  for (int j=0;j<DIM;++j)
	    b[0] += __Ftilde(j,i) * U[i+1][j+1] * N.m();   
	  
	
      }
    else if ((__domain>0))  // structure
      {
	// pressure: duale Laplace Fortsetzung
  	for (int j=0;j<DIM;++j)
  	  b[0] += __alpha_p * U[0][j+1] * N[j+1];  

      }
    else abort();

  }
  
  /*-----------------------------------------*/
  
  template<int DIM>
  void AleDual<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    /// not necessary for dual problem and be derived by using A^T
    std::cerr << "Matrix wird nicht gebraucht" << std::endl;
    abort();
    
  }


  ////////////////////////////////////////////////// LPS

  template<int DIM>
  void AleDual<DIM>::lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    __domain = __chi(v);
    
    double vm = 0;
    for (int d=0;d<DIM;++d) vm += U[d+1].m() * U[d+1].m();
    vm = sqrt(vm);
    //    vm = 0.0;
    
    __delta_lps_p =  __delta_lps_p0/__rho_f / (__nu_f/ (h*h) + vm / h );
    if (__domain>=0) __delta_lps_p = 0;

    __delta_lps_v =  __delta_lps_v0*__rho_f / (__nu_f/ (h*h) + vm / h );
    if (__domain>=0) __delta_lps_v = 0;
  }
  
    
  template<int DIM>
  void AleDual<DIM>::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
  {
    
    for (int j=0;j<DIM;++j)
      b[0] += __delta_lps_p * UP[0][j+1] * N[j+1];


    double conv_N = 0;
    for (int j=0;j<DIM;++j)
      conv_N += __delta_lps_v * U[j+1].m() * N[j+1];
    for (int d=0;d<DIM;++d)
      for (int j=0;j<DIM;++j)
	b[d+1] += U[j+1].m() * UP[d+1][j+1] * conv_N;
  
    
 }
  
  
  template<int DIM>
  void AleDual<DIM>::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
  {
    for (int j=0;j<DIM;++j)
      A(0,0) += __delta_lps_p * Mp[j+1] * Np[j+1];

    double conv_N = 0;
    for (int j=0;j<DIM;++j)
      conv_N += __delta_lps_v * U[j+1].m() * Np[j+1];
    double conv_M = 0;
    for (int j=0;j<DIM;++j)
      conv_M += U[j+1].m() * Mp[j+1];
    for (int d=0;d<DIM;++d)
      A(d+1,d+1) += conv_M * conv_N;

 }
  


template class AleDual<2>;
template class AleDual<3>;

}



/*-----------------------------------------*/
