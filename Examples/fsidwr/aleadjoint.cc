#include  "ale.h"
#include  "filescanner.h"
#include  "stopwatch.h"


/*-----------------------------------------*/

namespace Gascoigne
{
  

  template<int DIM>
  AleAdjoint<DIM>::AleAdjoint()
  {
    abort();
  }
  
     
  
  template<int DIM>
  AleAdjoint<DIM>::AleAdjoint(const ParamFile* pf) : Equation(), AleBase<DIM>()
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
    FileScanner FS(DFH,pf,"Equation");

    __chi.BasicInit(__solid_type);
  }


  template<int DIM>
  void AleAdjoint<DIM>::point(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    __domain = __chi(v);
    __h = h;
    __v = v;
    __alpha_p = __alpha_p0 * h * h / __nu_f;
    AleBase<DIM>::compute_transformation(U,__domain);
  }
  
  
  template<int DIM>
  void AleAdjoint<DIM>::Form(VectorIterator b, const FemFunction& U, 
			     const TestFunction& N) const
  {
    if (__domain>0)  // solid
      {
	// ------------------------------ Extension of pressure
	for (int i=0;i<DIM;++i)
	  b[0] += U[0][i+1] * N[i+1];

	// ---------------------------------------- TENSOR
	// -------------------- 2\mu F E
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int k=0;k<DIM;++k)	      
	      b[1+i] += 2.0 * __mu_s * 
		AleBase<DIM>::__F(i,k) * AleBase<DIM>::__E(k,j) * N[j+1];

	// -------------------- \lambda (trace E) F * \nabla\phi
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[1+i] +=__lambda_s * AleBase<DIM>::__trace_E * 
	      AleBase<DIM>::__F(i,j) * N[j+1];


	// -------------------- Extension of Velocity (zero)
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[1+DIM+i] += U[1+i][j+1] * N[j+1];
      }
    else             // fluid
      {
	for (int i=0;i<DIM;++i)
	  b[0] += __alpha_p *  U[0][i+1] * N[i+1];

	// -------------------- TENSOR - PRESSURE  (-p Ftilde^T, nabla\phi)
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[1+i] -= U[0].m() * AleBase<DIM>::__Ftilde(j,i) * N[j+1];

	// -------------------- TENSOR - VELOCITY (
	//   \rho\nu (nVF + nNF^T) Ftilde^T/J, \nabla \phi
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int k=0;k<DIM;++k)
	      b[1+i] += __rho_f * __nu_f/AleBase<DIM>::__J *
		(AleBase<DIM>::__nablaV_Ftilde(i,k)+AleBase<DIM>::__nablaV_Ftilde(k,i)) 
		* AleBase<DIM>::__Ftilde(j,k) 
		* N[j+1]; 
	

	// -------------------- Divergence div (JF^{-1} v)  = div (Ftilde v)
	for (int j=0;j<DIM;++j)
	  for (int k=0;k<DIM;++k)
	    b[0] += AleBase<DIM>::__Ftilde(j,k) * U[k+1][j+1] * N.m();

	// -------------------- Convection    (Ftilde V) \nabla V
	// -------------------- \nabla V Ftilde v
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[i+1] += __rho_f * AleBase<DIM>::__nablaV_Ftilde(i,j) * U[j+1].m() * N.m();
	

	// -------------------- Extension of the displacement
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    b[1+DIM+i] += U[1+DIM+i][j+1] * N[j+1];
      }
  }
  

  template class AleAdjoint<2>;
  template class AleAdjoint<3>;

}



/*-----------------------------------------*/
