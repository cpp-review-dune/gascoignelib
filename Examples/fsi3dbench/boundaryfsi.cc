#include  "boundaryfsi.h"
#include  "filescanner.h"

extern double __DT;
extern double __THETA;



using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

  ////////////////////////////////////////////////// 







  ////////////////////////////////////////////////// BOUNDARY


  


  void BoundaryFSI::Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
  {
   
    
		if (col==4)
		{	
			for (int i=0;i<3;++i)
				 b[i+1] += __n3d[i]*0*N.m();
		}
		
    if (col==2)
		{	
			for (int i=0;i<3;++i)
				 b[i+1] += __n3d[i]*0*N.m();
		}
  }
  

   void BoundaryFSI::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
   {
   }
	
  


   void BoundaryFSI::pointboundary(double h, const FemFunction& U, const Vertex3d& v, const Vertex3d& n) const
  {

 			 __n3d=n;
	}
  
  

  
  

  

}

/*-----------------------------------------*/
