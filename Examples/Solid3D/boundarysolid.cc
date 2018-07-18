#include  "boundarysolid.h"
#include  "filescanner.h"
#include "multiplexbound.xx"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

  template<int DIM>
  BoundarySolid<DIM>::BoundarySolid(const ParamFile* pf) 
  {
    DataFormatHandler DFH;
    DFH.insert("pressure" ,&__pressure);
    FileScanner FS(DFH,pf,"Equation");
  }

  template<int DIM>
  void BoundarySolid<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
  {

   // if (col==2)
   	if (col==1 ||col==2 ||col==3 ||col==4 )
		{	
			VECTOR g=J*__pressure*F.inverse().transpose()*__n;

			for (int i=0;i<DIM;++i)
				{
				 b[i] += g[i]*N.m();
				}
		}
		
  }

  template<int DIM>
   void BoundarySolid<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
   {

    //if (col==2)
   	if (col==1 ||col==2 ||col==3 ||col==4 )
		 {	
				VECTOR phi; multiplexbound_init_test<DIM>(phi,N);
				VECTOR psi; multiplexbound_init_test<DIM>(psi,M);


				for (int j=0;j<DIM;++j)
					{
						g_dU[j] = J*(F.inverse().block(0,j,DIM,1)*psi.transpose()).trace()*__pressure*F.inverse().transpose()*__n-J*__pressure*F.inverse().transpose()*psi*F.inverse().transpose().block(j,0,1,DIM)*__n;
					}
		
				for (int j=0;j<DIM;++j) // F Sigma_j \nabla phi
					{
						VECTOR X = g_dU[j]*N.m();
						for (int i=0;i<DIM;++i)
							{
								A(i,j) +=  X(i,0) ;
							}
					}
			}
   }
	
  

  template<int DIM>
   void BoundarySolid<DIM>::pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const
  {
 			 __n[0]=n[0];__n[1]=n[1]; if(DIM==3)__n[2]=n[2]; 
 			
 			multiplexbound_init_F<DIM>(F,U);
      J = F.determinant();
	}
  
    template<int DIM>
   void BoundarySolid<DIM>::point_M(int j, const FemFunction& U, const TestFunction& M) const
  {
  //wird nicht aufgerufen!!!
    /*VECTOR psi; multiplexbound_init_test<DIM>(psi,M);

			for (int j=0;j<DIM;++j)
				{
					g_dU[j] = J*(F.inverse().block(0,j,DIM,1)*psi.transpose()).trace()*__pressure*F.inverse().transpose()*__n-J*__pressure*F.inverse().transpose()*psi*F.inverse().transpose().block(j,0,1,DIM)*__n;
				}
	*/
  }
  



  
  template class BoundarySolid<2>;
  template class BoundarySolid<3>;

  

}

/*-----------------------------------------*/
