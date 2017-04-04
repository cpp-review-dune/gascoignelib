#include  "ale_slow.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
  template<int DIM>
  AleSlow<DIM>::AleSlow()
  {
    abort();
  }
  
  /*-----------------------------------------*/
  
  template<int DIM>
  AleSlow<DIM>::AleSlow(const ParamFile* pf) : Equation()
  {
    std::string __solid_type;
    
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
  void AleSlow<DIM>::compute_F(const FemFunction& U) const
  {
    __F.zero();
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	__F(i,j) = U[i+1+DIM][j+1];     // ok, ok wicki
    for (int i=0;i<DIM;++i) 
      __F(i,i) += 1.0;                  // ok, ok wicki
  }
  template<int DIM>
  void AleSlow<DIM>::compute_Ftilde() const
  {
    // Ftilde
    if (DIM==2)
      {
	__Ftilde(0,0) =  __F(1,1);
	__Ftilde(1,1) =  __F(0,0);
	__Ftilde(1,0) = -__F(1,0);
	__Ftilde(0,1) = -__F(0,1);          // ok , ok wicki
      }
    else
      {
	__Ftilde.zero();
	__Ftilde(0, 0)= __F(1,1)*__F(2,2)-__F(1,2)*__F(2,1);
	__Ftilde(0, 1)=-__F(0,1)*__F(2,2)+__F(0,2)*__F(2,1);
	__Ftilde(0, 2)= __F(0,1)*__F(1,2)-__F(0,2)*__F(1,1);
	__Ftilde(1, 0)=-__F(1,0)*__F(2,2)+__F(1,2)*__F(2,0);
	__Ftilde(1, 1)= __F(0,0)*__F(2,2)-__F(0,2)*__F(2,0);
	__Ftilde(1, 2)=-__F(0,0)*__F(1,2)+__F(0,2)*__F(1,0);
	__Ftilde(2, 0)= __F(1,0)*__F(2,1)-__F(1,1)*__F(2,0);
	__Ftilde(2, 1)=-__F(0,0)*__F(2,1)+__F(0,1)*__F(2,0);
	__Ftilde(2, 2)= __F(0,0)*__F(1,1)-__F(0,1)*__F(1,0);
      }
  }

  template<int DIM>
  void AleSlow<DIM>::compute_convection(const FemFunction& U) const
  {
    __conv_v=0;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	__conv_v[i] += U[j+1].m() * U[i+1][j+1]; // ok , ok wicki
  }


  template<int DIM>
  double AleSlow<DIM>::DU_trace_E(int k, const TestFunction& M) const
  {
    double res = 0.0;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	res += DU_F(i,j,k,M)*__F(i,j);
    return res;
  }

  template<int DIM>
  double AleSlow<DIM>::DU_F(int i,int j, int k, const TestFunction& M) const
  {
    if (i!=k) return 0.0;
    return M[j+1];
  }
  

  template<int DIM>
  double AleSlow<DIM>::DU_divergence_v(int k, const FemFunction& U, const TestFunction& M) const
  {
    double res = 0.0;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	res += DU_Ftilde(i,j,k,M) * U[j+1][i+1];
    return res;
  }
  template<int DIM>
  double AleSlow<DIM>::DV_divergence_v(int k, const TestFunction& M) const
  {
    double res = 0.0;
    for (int i=0;i<DIM;++i)
      res +=  __Ftilde(i,k) * M[i+1];
    return res;
  }
  

  template<int DIM>
  double AleSlow<DIM>::DU_J(int k, const FemFunction& U, const TestFunction& M) const
  {


    if (DIM==2)
      {
	if (k==0) return M.x() * (1.0+U[1+DIM+1].y()) - M.y() * U[1+DIM+1].x();
	if (k==1) return M.y() * (1.0+U[1+DIM+0].x()) - M.x() * U[1+DIM+0].y();
      }
    else if (DIM==3)
      {
#define U1 U[1+DIM+0]
#define U2 U[1+DIM+1]
#define U3 U[1+DIM+2]
#define UK U[1+DIM+k]
	
	double res = M[k+1];
	
	res += M[k+1] * ( U1.x()+U2.y()+U3.z() - UK[k+1]);
	res -= M.x() * U1[k+1] + M.y() * U2[k+1] + M.z() * U3[k+1] - M[k+1]*UK[k+1];
	
	if (k==0) 
	  {
	    res += M.x() * ( U2.y()*U3.z()-U2.z()*U3.y());
	    res += M.y() * (-U2.x()*U3.z()+U2.z()*U3.x());
	    res += M.z() * ( U2.x()*U3.y()-U2.y()*U3.x());
	  }
	if (k==1) 
	  {
	    res += M.x() * (-U1.y()*U3.z()+U1.z()*U3.y());
	    res += M.y() * ( U1.x()*U3.z()-U1.z()*U3.x());
	    res += M.z() * (-U1.x()*U3.y()+U1.y()*U3.x());
	  }
	if (k==2) 
	  {
	    res += M.x() * ( U1.y()*U2.z()-U1.z()*U2.y());
	    res += M.y() * (-U1.x()*U2.z()+U1.z()*U2.x());
	    res += M.z() * ( U1.x()*U2.y()-U1.y()*U2.x());
	  }
	//	"1+u3z+u2y+u2y*u3z-u2z*u3y+u1x+u1x*u3z+u1x*u2y+u1x*u2y*u3z-u1x*u2z*u3y-u2x*u1y-u2x*u1y*u3z+u2x*u1z*u3y+u3x*u1y*u2z-u3x*u1z-u3x*u1z*u2y"
	  
	return res;
	
#undef U1
#undef U2
#undef U3
#undef UK
      }
    
    abort();
    
  }
  


  template<int DIM>
  double AleSlow<DIM>::DU_nablaV_Ftilde(int i, int j, int k, const FemFunction& U, const TestFunction& M) const
  {
    if (DIM==3) return 0.0;
    
    double res = 0.0;
    for (int l=0;l<DIM;++l)
      res += U[i+1][l+1] * DU_Ftilde(l,j,k,M);
    return res;
  }
  template<int DIM>
  double AleSlow<DIM>::DV_nablaV_Ftilde(int i, int j, int k, const TestFunction& M) const 
  {
    if (i!=k) return 0.0;
    double res = 0.0;
    for (int l=0;l<DIM;++l)  res += M[l+1] * __Ftilde(l,j);
    return res;
  }

  // --------------------------------------------------

  template<int DIM>
  double AleSlow<DIM>::DV_conv_v(int i, int k, const FemFunction& U, const TestFunction& M) const
  {
    double res = M.m() * U[i+1][k+1];
    if (k==i)
      for (int j=0;j<DIM;++j) 
	res += U[j+1].m() * M[j+1];
    return res;
  }

  template<int DIM>
  double AleSlow<DIM>::DU_Ftilde(int i, int j, int k, const TestFunction& M) const
  {
    if (DIM==3) 
      {
	// das ist sehr ineffizient.... ich glaube, 
	//    --  wenn  j==k, dann ist das ergebnis immer null
	if ((i==0)&&(j==0)) return ( DU_F(1,1,k,M)*__F(2,2)-DU_F(1,2,k,M)*__F(2,1))+( __F(1,1)*DU_F(2,2,k,M)-__F(1,2)*DU_F(2,1,k,M));
	if ((i==0)&&(j==1)) return (-DU_F(0,1,k,M)*__F(2,2)+DU_F(0,2,k,M)*__F(2,1))+(-__F(0,1)*DU_F(2,2,k,M)+__F(0,2)*DU_F(2,1,k,M));
	if ((i==0)&&(j==2)) return ( DU_F(0,1,k,M)*__F(1,2)-DU_F(0,2,k,M)*__F(1,1))+( __F(0,1)*DU_F(1,2,k,M)-__F(0,2)*DU_F(1,1,k,M));
	if ((i==1)&&(j==0)) return (-DU_F(1,0,k,M)*__F(2,2)+DU_F(1,2,k,M)*__F(2,0))+(-__F(1,0)*DU_F(2,2,k,M)+__F(1,2)*DU_F(2,0,k,M));
	if ((i==1)&&(j==1)) return ( DU_F(0,0,k,M)*__F(2,2)-DU_F(0,2,k,M)*__F(2,0))+( __F(0,0)*DU_F(2,2,k,M)-__F(0,2)*DU_F(2,0,k,M));
	if ((i==1)&&(j==2)) return (-DU_F(0,0,k,M)*__F(1,2)+DU_F(0,2,k,M)*__F(1,0))+(-__F(0,0)*DU_F(1,2,k,M)+__F(0,2)*DU_F(1,0,k,M));
	if ((i==2)&&(j==0)) return ( DU_F(1,0,k,M)*__F(2,1)-DU_F(1,1,k,M)*__F(2,0))+( __F(1,0)*DU_F(2,1,k,M)-__F(1,1)*DU_F(2,0,k,M));
	if ((i==2)&&(j==1)) return (-DU_F(0,0,k,M)*__F(2,1)+DU_F(0,1,k,M)*__F(2,0))+(-__F(0,0)*DU_F(2,1,k,M)+__F(0,1)*DU_F(2,0,k,M));
	if ((i==2)&&(j==2)) return ( DU_F(0,0,k,M)*__F(1,1)-DU_F(0,1,k,M)*__F(1,0))+( __F(0,0)*DU_F(1,1,k,M)-__F(0,1)*DU_F(1,0,k,M));
	abort();
      }
    else if (DIM==2)
      {
	if ( (i==0)&&(j==0)&&(k==1) ) return  M.y();
	if ( (i==0)&&(j==1)&&(k==0) ) return -M.y();
	if ( (i==1)&&(j==0)&&(k==1) ) return -M.x();
	if ( (i==1)&&(j==1)&&(k==0) ) return  M.x();
	return 0.0;
      }    
    else abort();
  }

  /*-----------------------------------------*/

  template<int DIM>
  void AleSlow<DIM>::point_solid(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    compute_F(U);
    
    __trace_E = -0.5 * DIM;
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	__trace_E += 0.5 * __F(j,i)*__F(j,i);  // ok wicki
  }

  template<int DIM>
  void AleSlow<DIM>::point_fluid(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    compute_F(U);
    __J = __F.det();
    compute_Ftilde();
    compute_convection(U);
    
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
  }
  

  template<int DIM>
  void AleSlow<DIM>::point(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    __domain = __chi(v);


    if (__domain>0)        // solid
      point_solid(h,U,v);
    else if (__domain<0)   // fluid
      point_fluid(h,U,v);
    else abort();
    

    __alpha_v = __alpha_v0/h/h;


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
    __alpha_p =  __alpha_p0*h*h;



    double vm = 0;
    for (int d=0;d<DIM;++d) vm += U[d+1].m() * U[d+1].m();
    vm = sqrt(vm);
    vm = 0.0;
    
    __delta_lap_p =  __delta_lap_p0/__rho_f / (__nu_f/ (h*h));
    if (__domain>=0) __delta_lap_p = 0;
  }

  /*-----------------------------------------*/

  template<int DIM>
  void AleSlow<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    for (int i=0;i<2;++i)
      b[i+1+2] += U[i+1+2].x() * N.x() + U[i+1+2].y() * N.y();
    for (int i=0;i<2;++i)
      b[i+1+2] += - 0.01*N.m();
    
    
    if (__domain<0) // fluid
      {
	// stab
	for (int j=0;j<DIM;++j)
	  b[0] += __delta_lap_p * U[0][j+1] * N[j+1];
	
	// div
	b[0] += __divergence_v * N.m();  // ok , ok wicki

	// convection
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int k=0;k<DIM;++k)
	      b[i+1] += __rho_f * __Ftilde(j,k) * U[k+1].m() * U[i+1][j+1] * N.m();


	// tensor
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int k=0;k<DIM;++k)
	      b[i+1] += __rho_f * __nu_f * (__nablaV_Ftilde(i,k) 
					    + __nablaV_Ftilde(k,i)
					    )   // ok ??? , scheint ok zu sein
		* __Ftilde(j,k) * N[j+1]  / __J;

	// pressure
	for (int d=0;d<DIM;++d)
	  for (int j=0;j<DIM;++j)
	    b[d+1] -= U[0].m() * __Ftilde(j,d) * N[j+1];  // ok  , ok wicki
	
	
	// laplace u
	for (int d=0;d<DIM;++d)
	  {
	    for (int j=0;j<DIM;++j)
	      b[d+1+DIM] += __alpha_u * U[d+1+DIM][j+1] * N[j+1];   // ok wicki
	  }
	
	
      }
    else if ((__domain>0))  // structure
      {
	// pressure: Laplace Fortsetzung
  	for (int j=0;j<DIM;++j)
  	  b[0] += __alpha_p * U[0][j+1] * N[j+1];   // ok wicki

	// tensor

	// 1. mu  F F^T F
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int l=0;l<DIM;++l)
	      for (int m=0;m<DIM;++m)
		b[i+1] +=  __mu_s * __F(i,l) * __F(m,l) * __F(m,j) * N[j+1]; // ok (hinsichtlich der auftretenden Terme) wicki
	// 2. F (lambda spur - mu)
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[i+1] += (__trace_E*__lambda_s - __mu_s) * __F(i,j) * N[j+1]; // ok (hinsichtlich der auftretenden Terme) wicki
	  

	
// 	for (int d=0;d<DIM;++d)
// 	  for (int j=0;j<DIM;++j)
// 	    b[d+1] += __rho_s * __mu_s * U[d+1+DIM][j+1] * N[j+1];

	// v=0
	for (int d=0;d<DIM;++d)
	  b[d+1+DIM] -= __alpha_v * U[d+1].m() * N.m();   // minus Zeichen evtl. weg
      }
    else abort();
  }
  
  /*-----------------------------------------*/
  
  template<int DIM>
  void AleSlow<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    for (int i=0;i<2;++i)
      A(i+1+2,i+1+2) += M.x() * N.x() + M.y() * N.y();
    
    if (__domain<0) // fluid
      MatrixFluid(A,U,M,N);
    else if (__domain>0) // solid
      MatrixSolid(A,U,M,N);
    else abort();
  }

  
  template<int DIM>
  void AleSlow<DIM>::MatrixSolid(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    //	pressure
    for (int j=0;j<DIM;++j)
      A(0,0) += __alpha_p * M[j+1] * N[j+1];
    
    // tensor
    // 1. mu  F F^T F
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	for (int l=0;l<DIM;++l)
	  for (int m=0;m<DIM;++m)
	    for (int k=0;k<DIM;++k)
	      {
		A(i+1,k+1+DIM) += __mu_s * DU_F(i,l,k,M) * __F(m,l) * __F(m,j) * N[j+1];
		A(i+1,k+1+DIM) += __mu_s * __F(i,l) * DU_F(m,l,k,M) * __F(m,j) * N[j+1];
		A(i+1,k+1+DIM) += __mu_s * __F(i,l) * __F(m,l) * DU_F(m,j,k,M) * N[j+1];
	      }
    
    // 2. F (lambda spur - mu)
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	for (int k=0;k<DIM;++k)
	  A(i+1,k+1+DIM) += (__trace_E * __lambda_s - __mu_s) * DU_F(i,j,k,M) * N[j+1];
    // 2. F (lambda spur - mu)
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	for (int k=0;k<DIM;++k)
	  A(i+1,k+1+DIM) += DU_trace_E(k,M) * __lambda_s * __F(i,j) * N[j+1];
    
    // v=0
    for (int d=0;d<DIM;++d)
      A(d+1+DIM,d+1) -= __alpha_v * M.m() * N.m();
  }
  
  
  template<int DIM>
  void AleSlow<DIM>::MatrixFluid(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    double lap = 0;
    for (int i=0;i<DIM;++i) lap += M[i+1] * N[i+1];
    
    // stab
    for (int j=0;j<DIM;++j)
      A(0,0) += __delta_lap_p * lap;
    

    // div
    for (int k=0;k<DIM;++k)
      {
	A(0,k+1)     += DV_divergence_v(k,M)   * N.m();
	A(0,k+1+DIM) += DU_divergence_v(k,U,M) * N.m();
      }
    
    
    
    // convection
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	for (int k=0;k<DIM;++k)
	  {
	    A(i+1,i+1) += __rho_f * __Ftilde(j,k) * U[k+1].m() * M[j+1] * N.m();
	    A(i+1,k+1) += __rho_f * __Ftilde(j,k) * M.m() * U[i+1][j+1] * N.m();
	    for (int l=0;l<DIM;++l)
	      A(i+1,l+1+DIM) += __rho_f * DU_Ftilde(j,k,l,M) * U[k+1].m() * U[i+1][j+1] * N.m();
	  }
    
    
    // tensor
    for (int d=0;d<DIM;++d)
      for (int i=0;i<DIM;++i)
	for (int j=0;j<DIM;++j)
	  {
	    for (int k=0;k<DIM;++k)
	      {
		A(d+1,k+1) += __rho_f * __nu_f * (DV_nablaV_Ftilde(d,i,k,M)+ DV_nablaV_Ftilde(i,d,k,M)) 
		  * __Ftilde(j,i) * N[j+1]  / __J;
		
		A(d+1,k+1+DIM) += __rho_f * __nu_f * (DU_nablaV_Ftilde(d,i,k,U,M)+DU_nablaV_Ftilde(i,d,k,U,M)) 
		  * __Ftilde(j,i) * N[j+1]  / __J;
		
		A(d+1,k+1+DIM) += __rho_f * __nu_f * (__nablaV_Ftilde(d,i)+ __nablaV_Ftilde(i,d)) 
		  * DU_Ftilde(j,i,k,M) * N[j+1]  / __J;
		
		A(d+1,k+1+DIM) += __rho_f * __nu_f * (__nablaV_Ftilde(d,i)+ __nablaV_Ftilde(i,d)) 
		  * __Ftilde(j,i ) * N[j+1]  * (- DU_J(k,U,M))/__J/__J;
	      }
	  }
    
    // pressure
    for (int d=0;d<DIM;++d)
      for (int j=0;j<DIM;++j)
	{
	  A(d+1,0) -= M.m() * __Ftilde(j,d) * N[j+1];
	  for (int k=0;k<DIM;++k)
	    A(d+1,k+1+DIM) -= U[0].m() * DU_Ftilde(j,d,k,M) * N[j+1];
	}
    
    // laplace u
    for (int d=0;d<DIM;++d)
      for (int j=0;j<DIM;++j)
	A(d+1+DIM,d+1+DIM) += __alpha_u * lap;
  }


  ////////////////////////////////////////////////// LPS

  template<int DIM>
  void AleSlow<DIM>::lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const
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
  void AleSlow<DIM>::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
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
  void AleSlow<DIM>::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
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
  


template class AleSlow<2>;
template class AleSlow<3>;

}



/*-----------------------------------------*/
