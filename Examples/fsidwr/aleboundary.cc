#include  "aleboundary.h"
#include  "filescanner.h"


/*-----------------------------------------*/

namespace Gascoigne
{
  
  /*-----------------------------------------*/
  
  template<int DIM>
  AleBoundary<DIM>::AleBoundary(const ParamFile* pf) : BoundaryEquation(),
						       AleBase<DIM>(),
						       __p_left(0.0), __p_right(0.0)
  {
    DataFormatHandler DFH;
    DFH.insert("solid_type",&__solid_type);
    DFH.insert("p_left",&__p_left,0.0);
    DFH.insert("p_right",&__p_right,0.0);
    DFH.insert("nu_f",&__nu_f,0.0);
    DFH.insert("rho_f" ,&__rho_f);
    FileScanner FS(DFH,pf,"Equation");

    __chi.BasicInit(__solid_type);
  }
  
  /*-----------------------------------------*/
  template<int DIM>
  void AleBoundary<DIM>::pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const
  {
    __n = n;
    __h = h;
    __v = v;
    __domain = __chi(v);
    AleBase<DIM>::compute_transformation(U,__domain);
  }

  /*-----------------------------------------*/

  template<int DIM>
  void AleBoundary<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N,int col) const
  {
    //  - n \cdot \rho\nu nNF^T Ftilde^T/J, \phi
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	for (int k=0;k<DIM;++k)
	  b[1+i+DIM] -= __rho_f * __nu_f *  __n[j] * N.m() *
	    AleBase<DIM>::__nablaV_Ftilde(k,i) * 
	    AleBase<DIM>::__Ftilde(j,k) 
	    / AleBase<DIM>::__J;
    
    
    if (col==0) 
      b[1+DIM] += __n.x() * __p_left * N.m();
    // for (int i=0;i<DIM;++i)
    // 	for (int j=0;j<DIM;++j)
    // 	  b[i+1] += __n[j] * AleBase<DIM>::__Ftilde(j,i) * __p_left * N.m();
    
    if (col==1) 
      b[1+DIM] += __n.x() * __p_right * N.m();
    
      // for (int i=0;i<DIM;++i)
      // 	for (int j=0;j<DIM;++j)
      // 	  b[i+1] += __n[j] * AleBase<DIM>::__Ftilde(j,i) * __p_right * N.m();
  }
  
  /*-----------------------------------------*/
  
  template<int DIM>
  void AleBoundary<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
  {
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	for (int k=0;k<DIM;++k)
	  for (int d=0;d<DIM;++d)
	    {
	      A(1+i+DIM,d+1+DIM) -= __rho_f * __nu_f *  __n[j] * N.m() *
		AleBase<DIM>::DU_nablaV_Ftilde(k,i,d,U,M) * 
		AleBase<DIM>::__Ftilde(j,k) 
		/ AleBase<DIM>::__J;
	      
	      A(1+i+DIM,d+1) -= __rho_f * __nu_f *  __n[j] * N.m() *
		AleBase<DIM>::DV_nablaV_Ftilde(k,i,d,U,M) * 
		AleBase<DIM>::__Ftilde(j,k) 
		/ AleBase<DIM>::__J;
	      
	      A(1+i+DIM,d+1+DIM) -= __rho_f * __nu_f *  __n[j] * N.m() *
		AleBase<DIM>::__nablaV_Ftilde(k,i) * 
		AleBase<DIM>::DU_Ftilde(j,k,d,U,M) 
		/ AleBase<DIM>::__J;
	      
	      A(1+i+DIM,d+1+DIM) -= __rho_f * __nu_f *  __n[j] * N.m() *
		AleBase<DIM>::__nablaV_Ftilde(k,i) * 
		AleBase<DIM>::__Ftilde(j,k) 
		* AleBase<DIM>::DU_J(d,U,M) * (-1.0) /
		(AleBase<DIM>::__J*AleBase<DIM>::__J);
	    }
    
    
    // if (col==1) 
    //   for (int i=0;i<DIM;++i)
    // 	for (int j=0;j<DIM;++j)
    // 	  for (int d=0;d<DIM;++d)
    // 	    A(i+1+DIM,d+1+DIM) += __n[j] * AleBase<DIM>::DU_Ftilde(j,i,d,U,M) * __p_left * N.m();
    
    // if (col==0) 
    //   for (int i=0;i<DIM;++i)
    // 	for (int j=0;j<DIM;++j)
    // 	  for (int d=0;d<DIM;++d)
    // 	    A(i+1+DIM,d+1+DIM) += __n[j] * AleBase<DIM>::DU_Ftilde(j,i,d,U,M) * __p_right * N.m();
    
    

    // double EPS = 1.e-6;
    // DoubleVector F1(GetNcomp()),F2(GetNcomp());
    // F1.zero();
    // pointboundary(__h, U, __v,__n);
    // Form(F1.begin(),U,N,col);
    // for (int i=0;i<GetNcomp();++i)
    //   {
    // 	FemFunction UU = U;
    // 	UU[i].add(EPS,M);
    // 	pointboundary(__h, UU, __v,__n);
    // 	F2.zero();
    // 	Form(F2.begin(),UU,N,col);
    // 	for (int j=0;j<GetNcomp();++j)
    // 	  A(j,i) += (F2[j]-F1[j])/EPS;
    //   }
  }



  /*-----------------------------------------*/
  
  template<int DIM>
  AleBoundaryAdjoint<DIM>::AleBoundaryAdjoint(const ParamFile* pf) : BoundaryEquation()
  {
    DataFormatHandler DFH;
    DFH.insert("solid_type",&__solid_type);
    DFH.insert("nu_f",&__nu_f,0.0);
    FileScanner FS(DFH,pf,"Equation");
    
    __chi.BasicInit(__solid_type);
  }
  
  /*-----------------------------------------*/
  template<int DIM>
  void AleBoundaryAdjoint<DIM>::pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const 
  {
    __n = n;
    __domain = __chi(v);
  }

  /*-----------------------------------------*/

  template<int DIM>
  void AleBoundaryAdjoint<DIM>::Form(VectorIterator b, const FemFunction& Z, const TestFunction& N,int col) const 
  {
    abort();
    if ((col==0)||(col==8))
      {
	b[1] -= __nu_f * (__n.x() * N.x()) * Z[1].m();
	b[2] -= __nu_f * (__n.y() * N.x()) * Z[1].m();

	b[1] -= __nu_f * (__n.x() * N.y()) * Z[2].m();
	b[2] -= __nu_f * (__n.y() * N.y()) * Z[2].m();
      }    
  }
  
  /*-----------------------------------------*/
  
  template<int DIM>
  void AleBoundaryAdjoint<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
  {
    abort();
    if ((col==0)||(col==8))
      {
	A(1,1) -= __nu_f * (__n.x() * N.x()) * M.m();
	A(2,1) -= __nu_f * (__n.y() * N.x()) * M.m();

	A(1,2) -= __nu_f * (__n.x() * N.y()) * M.m();
	A(2,2) -= __nu_f * (__n.y() * N.y()) * M.m();
      }
  }
  


template class AleBoundary<2>;
template class AleBoundary<3>;

template class AleBoundaryAdjoint<2>;
template class AleBoundaryAdjoint<3>;

}



/*-----------------------------------------*/
