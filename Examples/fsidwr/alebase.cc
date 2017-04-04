#include "alebase.h"

using namespace std;


namespace Gascoigne
{

  template<int DIM>
  void AleBase<DIM>::compute_transformation(const FemFunction& U, int domain) const
  {
    __F.zero();
    // -------------------- Compute F
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	__F(i,j) = U[i+1+DIM][j+1];  
    for (int i=0;i<DIM;++i) 
      __F(i,i) += 1.0;               
    // -------------------- Compute J
    __J = __F.det();


    if (domain<0) // fluid
      {
	// -------------------- Compute F-Tilde
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

	// -------------------- J \nabla v F^{-1} = \nabla v Ftilde
	__nablaV_Ftilde.zero();
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int k=0;k<DIM;++k)
	      __nablaV_Ftilde(i,j) += U[i+1][k+1] * __Ftilde(k,j);
      }
    
    if (domain>0) // solid
      {
	// -------------------- E = 1/2 (F^t F - I)
	__E.zero();
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    for (int k=0;k<DIM;++k)
	      __E(i,j) += __F(k,i) * __F(k,j);
	for (int i=0;i<DIM;++i)
	  __E(i,i) -= 1.0;
	__E *= 0.5;

	// -------------------- trace(E)
	__trace_E = 0.0;
	for (int i=0;i<DIM;++i)
	  __trace_E += __E(i,i);
      }
  }
  
  //////////////////////////////////////////////////

  template<int DIM> AleBase<DIM>::AleBase()
  {
    __F.resize(DIM,DIM);
    __Ftilde.resize(DIM,DIM);
    __nablaV_Ftilde.resize(DIM,DIM);
    __E.resize(DIM,DIM);
  }

  // --------------------------------------------------

  template<int DIM> 
  double AleBase<DIM>::DU_Ftilde(int i, int j, int d,
				 const FemFunction& U,
				 const TestFunction& M) const
  {
    if (DIM==2)
      {
	if ( (i==0)&&(j==0)&&(d==1) ) return  M.y();
	if ( (i==0)&&(j==1)&&(d==0) ) return -M.y();
	if ( (i==1)&&(j==0)&&(d==1) ) return -M.x();
	if ( (i==1)&&(j==1)&&(d==0) ) return  M.x();
	return 0;
      }
    else if (DIM==3)
      {
	if (j==d) return 0.0;
	double res = 0.0;
	
	if (i==j) res += M[d+1];   ////// denn i=j=d ist schon ausgeschlossen.	
	if (i==d) res -= M[j+1];   ////// denn i=j=d ist schon ausgeschlossen.

	return res;
	

	// int dx = 1+DIM+0;
	// int dy = 1+DIM+1;
	// int dz = 1+DIM+2;
	
	// + quadratisce Terme

	/*

       [u1y u2z - u1z u2y     -u0y u2z + u0z u2y    u0y u1z - u0z u1y ]
       [                                                              ]
       [-u1x u2z + u1z u2x    u0x u2z - u2x u0z     -u1z u0x + u0z u1x]
       [                                                              ]
       [u1x u2y - u2x u1y     -u2y u0x + u0y u2x    u0x u1y - u1x u0y ]
	 */
	
      }
    else { abort(); }
  }

  
  template<int DIM>
  double AleBase<DIM>::DU_J(int d, const FemFunction& U, const TestFunction& M) const 
  {
    if (DIM==2)
      {
	if (d==0) return M.x() * (1.0+U[1+DIM+1].y()) - M.y() * U[1+DIM+1].x();
	if (d==1) return M.y() * (1.0+U[1+DIM+0].x()) - M.x() * U[1+DIM+0].y();
	abort();
      }
    else if (DIM==3) 
      {
	int d1 = (d+1)%3;
	int d2 = (d+2)%3;
	
	double res = M[d+1]* (U[1+DIM+d1][1+d1] + U[1+DIM+d2][1+d2]);
	res -= M[d1+1] * U[1+DIM+d1][d+1];
	res -= M[d2+1] * U[1+DIM+d2][d+1];

	// Plus Terme dritter Ordnung...
	res += M[d+1] * U[1+DIM+d1][1+d1] * U[1+DIM+d2][1+d2];
	res -= M[d+1] * U[1+DIM+d1][1+d2] * U[1+DIM+d2][1+d1];

	res += M[d1+1] * U[1+DIM+d1][1+d2] * U[1+DIM+d2][1+d];
	res -= M[d1+1] * U[1+DIM+d1][1+d] * U[1+DIM+d2][1+d2];

	res += M[d2+1] * U[1+DIM+d1][1+d] * U[1+DIM+d2][1+d1];
	res -= M[d2+1] * U[1+DIM+d1][1+d1] * U[1+DIM+d2][1+d];
	return res;
      }
  }
  
  template<int DIM>
  double AleBase<DIM>::DU_F(int i, int j, int d, const FemFunction& U, const TestFunction& M) const 
  {
    if (i==d) return M[j+1];
    return 0.0;
  }

  template<int DIM>
  double AleBase<DIM>::DU_E(int i, int j, int d, const FemFunction& U, const TestFunction& M) const 
  {
    double res = 0.0;


    // for (int k=0;k<DIM;++k)
    //   res  += DU_F(k,i,d,U,M) * __F(k,j) + __F(k,i) * DU_F(k,j,d,U,M);
    // Es ist DU_F \neq 0 nur wenn i==d (also k==d)

    res  += DU_F(d,i,d,U,M) * __F(d,j) + __F(d,i) * DU_F(d,j,d,U,M);
    return 0.5 * res;
  }
  
  template<int DIM>
  double AleBase<DIM>::DU_trace_E(int d, const FemFunction& U, const TestFunction& M) const 
  {
    double res = 0.0;
    // for (int i=0;i<DIM;++i)
    //   for (int k=0;k<DIM;++k)
    // 	res  += DU_F(k,i,d,U,M) * __F(k,i);

    // DU_F \neq 0 nur falls k==d
    for (int i=0;i<DIM;++i)
      res  += DU_F(d,i,d,U,M) * __F(d,i);
    return res;
  }


  template<int DIM>
  double AleBase<DIM>::DU_nablaV_Ftilde(int i, int j, int d, const FemFunction& U, const TestFunction& M) const  
  {
    double res = 0.0;

    ///
    for (int k=0;k<DIM;++k)
      res +=  U[i+1][k+1] * DU_Ftilde(k,j,d,U,M);
    for (int k=0;k<DIM;++k)
      res +=  U[i+1][k+1] * DU_Ftilde(k,j,d,U,M);
    
    return res;
  }

  template<int DIM>
  double AleBase<DIM>::DV_nablaV_Ftilde(int i, int j, int d, const FemFunction& U, const TestFunction& M) const  
  {
    double res = 0.0;

    if (i==d)
      for (int k=0;k<DIM;++k)
	res += M[k+1] * __Ftilde(k,j);

    return res;
  }

  

template class AleBase<2>;
template class AleBase<3>;

}

