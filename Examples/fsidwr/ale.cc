#include  "ale.h"
#include  "filescanner.h"
#include  "stopwatch.h"


/*-----------------------------------------*/

namespace Gascoigne
{
  

  template<int DIM>
  Ale<DIM>::Ale()
  {
    abort();
  }
  
     
  
  template<int DIM>
  Ale<DIM>::Ale(const ParamFile* pf) : Equation(), AleBase<DIM>()
  {
    std::string __solid_type;
    
    DataFormatHandler DFH;
    DFH.insert("solid_type" ,&__solid_type);
    DFH.insert("alpha_p" ,&__alpha_p0);
    DFH.insert("alpha_lps" ,&__alpha_lps0);
    DFH.insert("nu_f" ,&__nu_f);
    DFH.insert("rho_f" ,&__rho_f);
    DFH.insert("rho_s" ,&__rho_s);
    DFH.insert("mu_s" ,&__mu_s);
    DFH.insert("lambda_s" ,&__lambda_s);

    DFH.insert("lambda_m" ,&__lambda_m);
    DFH.insert("mu_m" ,&__mu_m);

    DFH.insert("extension", &__extension, "harmonic");
    
    FileScanner FS(DFH,pf,"Equation");

    __chi.BasicInit(__solid_type);
  }


  template<int DIM>
  void Ale<DIM>::point(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    __domain = __chi(v);
    __h = h;
    __v = v;
    __alpha_p = __alpha_p0 * h * h / __nu_f;
    AleBase<DIM>::compute_transformation(U,__domain);
  }
  

  //
  // Trial:             Test: 
  //         0:   P           0    xi
  //         1/2  V           1/2  wf
  //         3/4  U           3/4  phi

  
  template<int DIM>
  void Ale<DIM>::Form(VectorIterator b, const FemFunction& U, 
		      const TestFunction& N) const
  {
    if (__domain>0)  // solid
      {
	
	// ---------------------------------------- TENSOR
	// -------------------- 2\mu F E
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int k=0;k<DIM;++k)	      
	      b[1+i+DIM] += 2.0 * __mu_s * 
		AleBase<DIM>::__F(i,k) * AleBase<DIM>::__E(k,j) * N[j+1];

	// -------------------- \lambda (trace E) F * \nabla\phi
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[1+i+DIM] +=__lambda_s * AleBase<DIM>::__trace_E * 
	      AleBase<DIM>::__F(i,j) * N[j+1];

      }
    else             // fluid
      {
	for (int i=0;i<DIM;++i)
	  b[0] += __alpha_p *  U[0][i+1] * N[i+1];

	// -------------------- TENSOR - PRESSURE  (-p Ftilde^T, nabla\phi)
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[1+i+DIM] -= U[0].m() * Ftilde(j,i) * N[j+1];

	// -------------------- TENSOR - VELOCITY (
	//   \rho\nu (nVF + nNF^T) Ftilde^T/J, \nabla \phi
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int k=0;k<DIM;++k)
	      b[1+i+DIM] += __rho_f * __nu_f/J() *
		(nablaV_Ftilde(i,k)+nablaV_Ftilde(k,i)) 
		* Ftilde(j,k) 
		* N[j+1]; 
	

	// -------------------- Divergence div (JF^{-1} v)  = div (Ftilde v)
	for (int j=0;j<DIM;++j)
	  for (int k=0;k<DIM;++k)
	    b[0] += Ftilde(j,k) * U[k+1][j+1] * N.m();

	// -------------------- Convection    (Ftilde V) \nabla V
	// -------------------- \nabla V Ftilde v
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[i+1+DIM] += __rho_f * AleBase<DIM>::__nablaV_Ftilde(i,j) * U[j+1].m() * N.m();
	

	// -------------------- Extension of the displacement
	if (__extension=="harmomic")
	  {
	    for (int i=0;i<DIM;++i)
	      for (int j=0;j<DIM;++j)
		b[1+i] += U[1+DIM+i][j+1] * N[j+1];
	  }
	else if (__extension == "elastic")
	  {
	    for (int i=0;i<DIM;++i)
	      for (int j=0;j<DIM;++j)
		{
		  b[1+i] += __mu_m * (U[1+DIM+i][j+1]+U[1+DIM+j][i+1])* N[j+1];
		  b[1+i] += __lambda_m * U[1+DIM+j][j+1] * N[i+1];
		}
	  }
	
      }
  }




  template<int DIM>
  void Ale<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    
    if (__domain>0)  // solid
      {
    	// ------------------------------ Extension of pressure
	// for (int i=0;i<DIM;++i)
	//   A(0,0) += M[i+1] * N[i+1];

    	// ---------------------------------------- TENSOR
    	// -------------------- 2\mu F E
    	for (int i=0;i<DIM;++i)
    	  for (int j=0;j<DIM;++j)
    	    for (int k=0;k<DIM;++k)	      
    	      for (int d=0;d<DIM;++d)	      
    		A(1+i+DIM,1+DIM+d) += 2.0 * __mu_s * N[j+1] *
    		  (AleBase<DIM>::DU_F(i,k,d,U,M) * AleBase<DIM>::__E(k,j)+
    		   AleBase<DIM>::__F(i,k) * AleBase<DIM>::DU_E(k,j,d,U,M));
	
		
    	// -------------------- \lambda (trace E) F * \nabla\phi
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int d=0;d<DIM;++d)
	      {
		A(1+i+DIM,1+d+DIM) +=__lambda_s * N[j+1] * 
		  (AleBase<DIM>::DU_trace_E(d,U,M) * AleBase<DIM>::__F(i,j) +
		   AleBase<DIM>::__trace_E * AleBase<DIM>::DU_F(i,j,d,U,M));
	      }
	
	

    	// // -------------------- Extension of Velocity (zero)
    	// for (int i=0;i<DIM;++i)
    	//   for (int j=0;j<DIM;++j)
    	//     A(1+DIM+i,1+i) += M[j+1] * N[j+1];
      }
    else             // fluid
      {
    	for (int i=0;i<DIM;++i)
    	  A(0,0) += __alpha_p *  M[i+1] * N[i+1];

    	// -------------------- TENSOR - PRESSURE  (-p Ftilde^T, nabla\phi)
    	for (int i=0;i<DIM;++i)
    	  for (int j=0;j<DIM;++j)
    	    {
    	      A(1+i+DIM,0) -= M.m() * Ftilde(j,i) * N[j+1];
    	      for (int d=0;d<DIM;++d)
    		A(1+i+DIM,1+DIM+d) -= U[0].m() * AleBase<DIM>::DU_Ftilde(j,i,d,U,M) * N[j+1];
    	    }

    	// -------------------- TENSOR - VELOCITY (
    	//   \rho\nu (nVF + nNF^T) Ftilde^T/J, \nabla \phi

    	for (int i=0;i<DIM;++i)
    	  for (int j=0;j<DIM;++j)
    	    for (int k=0;k<DIM;++k)
    	      {
    		// for (int d=0;d<DIM;++d)
    		//   A(1+i,1+d) += __rho_f * __nu_f/J() *
    		//     (AleBase<DIM>::DV_nablaV_Ftilde(i,k,d,U,M)+
    		//      AleBase<DIM>::DV_nablaV_Ftilde(k,i,d,U,M)) 
    		//     * Ftilde(j,k) 
    		//   * N[j+1]; 
		//    ##########     ---->>>  DV_nablaV_Ftilde \neq 0  i==d (bzw k==d)
		A(1+i+DIM,1+i) += __rho_f * __nu_f/J() *
		  (AleBase<DIM>::DV_nablaV_Ftilde(i,k,i,U,M))
		  * Ftilde(j,k) 
    		  * N[j+1]; 
		A(1+i+DIM,1+k) += __rho_f * __nu_f/J() *
		  (AleBase<DIM>::DV_nablaV_Ftilde(k,i,k,U,M)) 
		  * Ftilde(j,k) 
    		  * N[j+1]; 
		

    		for (int d=0;d<DIM;++d)
    		  A(1+i+DIM,1+d+DIM) += __rho_f * __nu_f/J() *
    		    (AleBase<DIM>::DU_nablaV_Ftilde(i,k,d,U,M)+
    		     AleBase<DIM>::DU_nablaV_Ftilde(k,i,d,U,M)) 
    		    * Ftilde(j,k) 
    		  * N[j+1]; 

    		for (int d=0;d<DIM;++d)
    		  A(1+i+DIM,1+d+DIM) += __rho_f * __nu_f/J() *
    		    (nablaV_Ftilde(i,k)+nablaV_Ftilde(k,i))
    		    * AleBase<DIM>::DU_Ftilde(j,k,d,U,M) 
    		    * N[j+1]; 

    		for (int d=0;d<DIM;++d)
    		  A(1+i+DIM,d+1+DIM) += __rho_f * __nu_f * AleBase<DIM>::DU_J(d,U,M) *
    		    (nablaV_Ftilde(i,k)+nablaV_Ftilde(k,i))
    		    * Ftilde(j,k) 
    		    * N[j+1] * (-1.0)/(J()*J());
    	      }
	


	
	
    	// -------------------- Divergence div (JF^{-1} v)  = div (Ftilde v)
    	for (int j=0;j<DIM;++j)
    	  for (int k=0;k<DIM;++k)
    	    {
    	      A(0,k+1) += Ftilde(j,k) * M[j+1] * N.m();
    	      for (int d=0;d<DIM;++d)
    		A(0,d+1+DIM) += AleBase<DIM>::DU_Ftilde(j,k,d,U,M) * U[k+1][j+1] * N.m(); 
    	    }
	
    	// -------------------- Convection    (Ftilde V) \nabla V
    	// -------------------- \nabla V Ftilde v
    	for (int i=0;i<DIM;++i)
    	  for (int j=0;j<DIM;++j)
    	    {
    	      for (int d=0;d<DIM;++d)
    		A(i+1+DIM,d+1+DIM) += __rho_f *  N.m() *
    		  AleBase<DIM>::DU_nablaV_Ftilde(i,j,d,U,M) * 
    		  U[j+1].m() ;
    	      for (int d=0;d<DIM;++d)
    		A(i+1+DIM,d+1) += __rho_f *  N.m() *
    		  AleBase<DIM>::DV_nablaV_Ftilde(i,j,d,U,M) * 
    		  U[j+1].m() ;
    	      A(i+1+DIM,j+1) += __rho_f *  N.m() * 
    		AleBase<DIM>::__nablaV_Ftilde(i,j) * 
    		M.m() ;
    	    }
	

	

    	// -------------------- Extension of the displacement
	if (__extension == "harmonic")
	  {
	    for (int i=0;i<DIM;++i)
	      for (int j=0;j<DIM;++j)
		A(1+i,1+DIM+i) += M[j+1] * N[j+1];
	  }
	else if (__extension == "elastic")
	  {
	    for (int i=0;i<DIM;++i)
	      for (int j=0;j<DIM;++j)
		{
		  A(1+i,1+DIM+i) += __mu_m * M[j+1] * N[j+1];
		  A(1+i,1+DIM+j) += __mu_m * M[i+1] * N[j+1];
		  A(1+i,1+DIM+j) += __lambda_m * M[j+1] * N[i+1];
		}
	  }
	else abort();

      }

    return;
    


    // double EPS = 1.e-6;
    // DoubleVector F1(GetNcomp()),F2(GetNcomp());
    // F1.zero();
    // point(__h, U, __v);
    // Form(F1.begin(),U,N);
    // for (int i=0;i<GetNcomp();++i)
    //   {
    // 	FemFunction UU = U;
    // 	UU[i].add(EPS,M);
    // 	point(__h, UU, __v);
    // 	F2.zero();
    // 	Form(F2.begin(),UU,N);
    // 	  for (int j=0;j<GetNcomp();++j)
    // 	    A(j,i) += (F2[j]-F1[j])/EPS;
    //   }
  }
  


  












  template class Ale<2>;
  template class Ale<3>;


}



/*-----------------------------------------*/
