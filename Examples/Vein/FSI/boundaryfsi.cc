#include  "boundaryfsi.h"
#include  "filescanner.h"

extern double __DT;
extern double __THETA;
extern double __TIME;

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

  ////////////////////////////////////////////////// 







  ////////////////////////////////////////////////// BOUNDARY
  template<int DIM>
  BoundaryFSI<DIM>::BoundaryFSI(const ParamFile* pf) 
  {
    DataFormatHandler DFH;
    DFH.insert("nu_f",&__nu_f,0.0);
    DFH.insert("rho_f" ,&__rho_f);
    FileScanner FS(DFH,pf,"Equation");
    p_2=2.266*1.0e4;
    p_4=2.286*1.0e4;
    	cout<<"%%%%%%%%%%Fluid_Stat%%%%%%%%%%"<<endl;
		cout<<"Boundary 2/4 -- do-nothing with p="<<p_2<<"g/cm/s^2="<<p_2/1333.22<<"mmHg and p="<<p_4<<"g/cm/s^2="<<p_4/1333.22<<"mmHg"<<endl;
		cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
  }

  template<int DIM>
  void BoundaryFSI<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
  {
  	//______________________________________________________________________________
  	if(DIM==3)
 			{
	 			NU << 
		    (*DEF)[0].x(), (*DEF)[0].y(),(*DEF)[0].z(),
		    (*DEF)[1].x(), (*DEF)[1].y(),(*DEF)[1].z(),
		    (*DEF)[2].x(), (*DEF)[2].y(),(*DEF)[2].z(),
		  	NV << 
		    U[1].x(), U[1].y(),U[1].z(),
		    U[2].x(), U[2].y(),U[2].z(),
		    U[3].x(), U[3].y(),U[3].z();
		  }
    //Berechnung einiger Werte
    // Deformationsgradient
    F= MATRIX::Identity()+NU;
    J = F.determinant();
    //______________________________________________________________________________
    if(DIM==3)
 			{
				NU_OLD << 
				  (*DEFOLD)[0].x(), (*DEFOLD)[0].y(),(*DEFOLD)[0].z(),
				  (*DEFOLD)[1].x(), (*DEFOLD)[1].y(),(*DEFOLD)[1].z(),
				  (*DEFOLD)[2].x(), (*DEFOLD)[2].y(),(*DEFOLD)[2].z(),
				NV_OLD << 
				  (*OLD)[1].x(), (*OLD)[1].y(),(*OLD)[1].z(),
				  (*OLD)[2].x(), (*OLD)[2].y(),(*OLD)[2].z(),
				  (*OLD)[3].x(), (*OLD)[3].y(),(*OLD)[3].z();
			}
    //Berechnung einiger Werte
    // Deformationsgradient
    F_OLD= MATRIX::Identity()+NU_OLD;
    J_OLD = F_OLD.determinant();
		//______________________________________________________________________________

    g 	  = -__rho_f * __nu_f * J*F.inverse().transpose()*NV.transpose()*F.inverse().transpose()*__n;
    g_OLD = -__rho_f * __nu_f * J_OLD*F_OLD.inverse().transpose()*NV_OLD.transpose()*F_OLD.inverse().transpose()*__n;
   

    double fact_time = 1.0;

    //if (__TIME<1.0) fact_time *= 0.5*(1.0-cos(M_PI*__TIME));

		if (col==4)
		{	
			for (int i=0;i<DIM;++i)
				{
				 b[i+1] += __THETA*g(i)*N.m()+(1.0-__THETA)*g_OLD(i)*N.m() +__n[i]*fact_time*(p_2)*N.m();
				 if((U[1].m()*__n[0]+U[2].m()*__n[1]+U[3].m()*__n[2])<0)
				 		b[i+1] -= 0.5*__THETA*(U[1].m()*__n[0]+U[2].m()*__n[1]+U[3].m()*__n[2])*U[1+i].m()*N.m();
				 if(((*OLD)[1].m()*__n[0]+(*OLD)[2].m()*__n[1]+(*OLD)[3].m()*__n[2])<0)
				 		b[i+1] -= 0.5*(1.0-__THETA)*((*OLD)[1].m()*__n[0]+(*OLD)[2].m()*__n[1]+(*OLD)[3].m()*__n[2])*(*OLD)[1+i].m()*N.m();
				 }
		}
		
    if (col==2)
		{	
			for (int i=0;i<DIM;++i)
				{
					 b[i+1] += __THETA*g(i)*N.m()+(1.0-__THETA)*g_OLD(i)*N.m()+__n[i]*fact_time*p_4*N.m();
				 if((U[1].m()*__n[0]+U[2].m()*__n[1]+U[3].m()*__n[2])<0)
				 		b[i+1] -= 0.5*__THETA*(U[1].m()*__n[0]+U[2].m()*__n[1]+U[3].m()*__n[2])*U[1+i].m()*N.m();
				 if(((*OLD)[1].m()*__n[0]+(*OLD)[2].m()*__n[1]+(*OLD)[3].m()*__n[2])<0)
				 		b[i+1] -= 0.5*(1.0-__THETA)*((*OLD)[1].m()*__n[0]+(*OLD)[2].m()*__n[1]+(*OLD)[3].m()*__n[2])*(*OLD)[1+i].m()*N.m();
				 }
		}
		
  }
  
  template<int DIM>
   void BoundaryFSI<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
   {
				 
			//_______________________________________________________________
			if(DIM==3)	
			{
				NU << 
				  (*DEF)[0].x(), (*DEF)[0].y(),(*DEF)[0].z(),
				  (*DEF)[1].x(), (*DEF)[1].y(),(*DEF)[1].z(),
				  (*DEF)[2].x(), (*DEF)[2].y(),(*DEF)[2].z(),
				NV << 
				  U[1].x(), U[1].y(),U[1].z(),
				  U[2].x(), U[2].y(),U[2].z(),
				  U[3].x(), U[3].y(),U[3].z();
			}
			//________________________________________________________________


					for (int j=0;j<DIM;++j)
							 {
							   NPHI=MATRIX::Zero();
							   PHI=VECTOR::Zero();
					 				for (int i=0;i<DIM;++i)
							   { 
						 			 NPHI(j,i)=M[i+1];
						 			 PHI(j)=M.m();
							   }
								//_________________________________________   
									g= - __rho_f * __nu_f*J*F.inverse().transpose()*NPHI.transpose()*F.inverse().transpose()*__n;
								         
								//___________________________________________
								for(int i=0;i<DIM;++i)
				 	         {
				 						A(1+i,j+1) +=__THETA*g(i)*N.m();					
				 						
				 						if((U[1].m()*__n[0]+U[2].m()*__n[1]+U[3].m()*__n[2])<0)
				 								{
				 									A(1+i,j+1)  -= 0.5*__THETA*(U[1].m()*__n[0]+U[2].m()*__n[1]+U[3].m()*__n[2])*PHI[i]*N.m();
				 								 	A(1+i,j+1) 	-= 0.5*__THETA*(PHI[0]*__n[0]+PHI[1]*__n[1]+PHI[2]*__n[2])*U[1+i].m()*N.m();      	 
				 								}
								   }
							  
				}
   }
	
  

  template<int DIM>
   void BoundaryFSI<DIM>::pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const
  {

 			 __n[0]=n[0];__n[1]=n[1]; if(DIM==3)__n[2]=n[2]; 
	}
  
  

  
  template class BoundaryFSI<2>;
  template class BoundaryFSI<3>;

  

}

/*-----------------------------------------*/
