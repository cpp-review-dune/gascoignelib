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
  void BoundaryFSI<DIM>::Form(VectorIterator b, const FemFunction& U_Dummy, const TestFunction& N, int col) const
  {
  	//______________________________________________________________________________
  	if(DIM==3)
 			{
	 			NU << 
		    (*U_Vec)[DIM+1+0].x(), (*U_Vec)[DIM+1+0].y(),(*U_Vec)[DIM+1+0].z(),
		    (*U_Vec)[DIM+1+1].x(), (*U_Vec)[DIM+1+1].y(),(*U_Vec)[DIM+1+1].z(),
		    (*U_Vec)[DIM+1+2].x(), (*U_Vec)[DIM+1+2].y(),(*U_Vec)[DIM+1+2].z();
		  	NV << 
		    (*U_Vec)[1+0].x(), (*U_Vec)[1+0].y(),(*U_Vec)[1+0].z(),
		    (*U_Vec)[1+1].x(), (*U_Vec)[1+1].y(),(*U_Vec)[1+1].z(),
		    (*U_Vec)[1+2].x(), (*U_Vec)[1+2].y(),(*U_Vec)[1+2].z();
		  }
    //Berechnung einiger Werte
    // Deformationsgradient
    F= MATRIX::Identity()+NU;
    J = F.determinant();
    //______________________________________________________________________________
    if(DIM==3)
 			{
	 			NU_OLD << 
				(*UOLD_Vec)[DIM+1+0].x(), (*UOLD_Vec)[DIM+1+0].y(),(*UOLD_Vec)[DIM+1+0].z(),
				(*UOLD_Vec)[DIM+1+1].x(), (*UOLD_Vec)[DIM+1+1].y(),(*UOLD_Vec)[DIM+1+1].z(),
				(*UOLD_Vec)[DIM+1+2].x(), (*UOLD_Vec)[DIM+1+2].y(),(*UOLD_Vec)[DIM+1+2].z();
				NV_OLD << 
				(*UOLD_Vec)[1+0].x(), (*UOLD_Vec)[1+0].y(),(*UOLD_Vec)[1+0].z(),
				(*UOLD_Vec)[1+1].x(), (*UOLD_Vec)[1+1].y(),(*UOLD_Vec)[1+1].z(),
				(*UOLD_Vec)[1+2].x(), (*UOLD_Vec)[1+2].y(),(*UOLD_Vec)[1+2].z();

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
				 if(((*U_Vec)[1+0].m()*__n[0]+(*U_Vec)[1+1].m()*__n[1]+(*U_Vec)[1+2].m()*__n[2])<0)
				 		b[i+1] -= 0.5*__THETA*((*U_Vec)[1+0].m()*__n[0]+(*U_Vec)[1+1].m()*__n[1]+(*U_Vec)[1+2].m()*__n[2])*(*U_Vec)[1+i].m()*N.m();
				 if(((*UOLD_Vec)[1+0].m()*__n[0]+(*UOLD_Vec)[1+1].m()*__n[1]+(*UOLD_Vec)[1+2].m()*__n[2])<0)
				 		b[i+1] -= 0.5*(1.0-__THETA)*((*UOLD_Vec)[1+0].m()*__n[0]+(*UOLD_Vec)[1+1].m()*__n[1]+(*UOLD_Vec)[1+2].m()*__n[2])*(*UOLD_Vec)[1+i].m()*N.m();
				 }
		}
		
    if (col==2)
		{	
			for (int i=0;i<DIM;++i)
				{
					 b[i+1] += __THETA*g(i)*N.m()+(1.0-__THETA)*g_OLD(i)*N.m()+__n[i]*fact_time*p_4*N.m();
				 if(((*U_Vec)[1+0].m()*__n[0]+(*U_Vec)[1+1].m()*__n[1]+(*U_Vec)[1+2].m()*__n[2])<0)
				 		b[i+1] -= 0.5*__THETA*((*U_Vec)[1+0].m()*__n[0]+(*U_Vec)[1+1].m()*__n[1]+(*U_Vec)[1+2].m()*__n[2])*(*U_Vec)[1+i].m()*N.m();
				 if(((*UOLD_Vec)[1+0].m()*__n[0]+(*UOLD_Vec)[1+1].m()*__n[1]+(*UOLD_Vec)[1+2].m()*__n[2])<0)
				 		b[i+1] -= 0.5*(1.0-__THETA)*((*UOLD_Vec)[1+0].m()*__n[0]+(*UOLD_Vec)[1+1].m()*__n[1]+(*UOLD_Vec)[1+2].m()*__n[2])*(*UOLD_Vec)[1+i].m()*N.m();
				 }
		}
		
  }
  
  template<int DIM>
   void BoundaryFSI<DIM>::Matrix(EntryMatrix& A, const FemFunction& U_Dummy, const TestFunction& M, const TestFunction& N, int col) const
   {
				 
			//_______________________________________________________________
			if(DIM==3)	
			{
	 			NU << 
				(*U_Vec)[DIM+1+0].x(), (*U_Vec)[DIM+1+0].y(),(*U_Vec)[DIM+1+0].z(),
				(*U_Vec)[DIM+1+1].x(), (*U_Vec)[DIM+1+1].y(),(*U_Vec)[DIM+1+1].z(),
				(*U_Vec)[DIM+1+2].x(), (*U_Vec)[DIM+1+2].y(),(*U_Vec)[DIM+1+2].z(),
			  	NV << 
				(*U_Vec)[1+0].x(), (*U_Vec)[1+0].y(),(*U_Vec)[1+0].z(),
				(*U_Vec)[1+1].x(), (*U_Vec)[1+1].y(),(*U_Vec)[1+1].z(),
				(*U_Vec)[1+2].x(), (*U_Vec)[1+2].y(),(*U_Vec)[1+2].z();
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
				 						
				 						 if(((*U_Vec)[1+0].m()*__n[0]+(*U_Vec)[1+1].m()*__n[1]+(*U_Vec)[1+2].m()*__n[2])<0)
				 								{
				 									A(1+i,j+1)  -= 0.5*__THETA*((*U_Vec)[1+0].m()*__n[0]+(*U_Vec)[1+1].m()*__n[1]+(*U_Vec)[1+2].m()*__n[2])*PHI[i]*N.m();
				 								 	A(1+i,j+1) 	-= 0.5*__THETA*(PHI[0]*__n[0]+PHI[1]*__n[1]+PHI[2]*__n[2])*(*U_Vec)[1+i].m()*N.m();      	 
				 								}
								   }
							  
				}
   }
	
  

  template<int DIM>
   void BoundaryFSI<DIM>::pointboundary(double h, const FemFunction& U_Dummy, const Vertex<DIM>& v, const Vertex<DIM>& n) const
  {

 			 __n[0]=n[0];__n[1]=n[1]; if(DIM==3)__n[2]=n[2]; 
	}
  
  

  
  template class BoundaryFSI<2>;
  template class BoundaryFSI<3>;

  

}

/*-----------------------------------------*/
