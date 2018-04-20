#include  "solideuler.h"
#include  "filescanner.h"



using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  template<int DIM>
  SolidEuler<DIM>::SolidEuler(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("rho_s" ,    &rho_s , 0.0);
    DFH.insert("lambda_s" ,    &lambda_s , 0.0);
    DFH.insert("mu_s" ,    &mu_s , 0.0);
    DFH.insert("mat_law" ,    &mat_law , "STVK");
    DFH.insert("rho_f" ,    &rho_f , 0.0);
    DFH.insert("nu_f" ,    &nu_f , 0.0);
    DFH.insert("lps" ,    &lps0 , 0.0);
       
    FileScanner FS(DFH, pf, "Equation");
    assert(rho_s>0);
    assert(lambda_s>0);
    assert(mu_s>0);
  
		kapa_s=lambda_s+2.0/3.0*mu_s;
    cout << "%%%%%%%%%% Problem %%%%%%%%%%" << endl
	  << "  rho_s / mu_s / lambda_s: " << rho_s << " / " << mu_s << " / " << lambda_s  << endl;
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;
    
  }

  ////////////////////////////////////////////////// 
 /*--------------------------------------------------------------------*/
 template<int DIM>
  void SolidEuler<DIM>::point_cell(int material) const 
{
		if(material==1) domain=1;
		if(material==2) domain=-1;
}
/*--------------------------------------------------------------------*/
#include "multiplex_solid_euler.xx"

 template<int DIM>
  void SolidEuler<DIM>::point(double h, const FemFunction& U, const Vertex<DIM>& v) const
  {
    //__h = h;
    //__v = v;

    if (domain<0)
      {//fluid
      	 multiplex_solid_euler_init_NV<DIM>(NV,*VEL);
		 multiplex_solid_euler_init_V<DIM>(V,*VEL);

		// TENSOR
		SIGMAf = rho_f * nu_f * ( NV + NV.transpose() );	
		
      }
    else
    {
		multiplex_solid_euler_init_f<DIM>(__f,U);

		__j = __f.determinant();
		
		if(mat_law=="STVK")
		{
			__e = 0.5*(__f.inverse().transpose()*__f.inverse()-MATRIX::Identity());
			SIGMAs = (2.0 * mu_s * __e + lambda_s * __e.trace() * MATRIX::Identity());
			}
			else if(mat_law=="artery")
		  {
		  __c=__f.inverse().transpose()*__f.inverse();
		  SIGMAs = mu_s*pow(__j,2.0/3.0)*(MATRIX::Identity()-1.0/3.0*__c.trace()*__c.inverse())+0.5*kapa_s*(pow(__j,-2.0)-1.0)*__c.inverse();
		  }
			else abort();
			sigmas=__j*__f.inverse()*SIGMAs*__f.inverse().transpose();
	}	
		
  }
  
  
  template<int DIM>
  void SolidEuler<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
      //phi =nabla N;
      VECTOR phi; multiplex_solid_euler_init_test<DIM>(phi,N);
      
      if (domain<0) // fluid
      {
      
		// divergence
		//b[0] += rho_f * NV.trace() * N.m();

		// tensor
		VECTOR X =  SIGMAf  * phi;
		for (int i=0;i<DIM;++i)
		  b[i+1] += X(i,0);

		// convection
		X =  rho_f * NV *V*N.m();
		for (int i=0;i<DIM;++i)
		  b[i+1] += X(i,0);


		// pressure
		X = -(*VEL)[0].m() * MATRIX::Identity()*phi;

		for (int i=0;i<DIM;++i)
		  b[i+1] += X(i,0);
		
      }
      if (domain>0) 
      {
		// FULL tensor F Sigma
		for (int i=0;i<DIM;++i) b[i+1] += (sigmas*phi)(i,0);
      }
  }
  

  template<int DIM>
  void SolidEuler<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {    
    if (domain>0) // solid
      {
		VECTOR phi; multiplex_solid_euler_init_test<DIM>(phi,N);
		VECTOR psi; multiplex_solid_euler_init_test<DIM>(psi,M);
    
		// F Sigma
		// wrt F, 	
		for (int j=0;j<DIM;++j) // F Sigma_j \nabla phi
			{
			  VECTOR X = sigmas_dU[j]*phi;
			  for (int i=0;i<DIM;++i)
			    {
			    	A(i+1,j+1) +=  X(i,0) ;
			    }
	  	}

	 }
    
  }



 template<int DIM>
 void SolidEuler<DIM>::point_M(int j, const FemFunction& U, const TestFunction& M) const
  {
      if (domain>0) // solid
      {
    	VECTOR psi; multiplex_solid_euler_init_test<DIM>(psi,M);
		for (int jj=0;jj<DIM;++jj)
			{
			 if(mat_law=="STVK")
		  	{
					MATRIX ej=  0.5 * (__f.inverse().transpose()*psi*__f.inverse().transpose().block(jj,0,1,DIM)*__f.inverse() +__f.inverse().transpose()*__f.inverse().block(0,jj,DIM,1)*psi.transpose()*__f.inverse());
					SIGMAs_dU= (2.0*mu_s*ej + lambda_s * ej.trace()*MATRIX::Identity());
				 }
	  	 else if(mat_law=="artery")
				{ 	
				    MATRIX __c_dU=  __f.inverse().transpose()*psi*__f.inverse().transpose().block(jj,0,1,DIM)*__f.inverse() +__f.inverse().transpose()*__f.inverse().block(0,jj,DIM,1)*psi.transpose()*__f.inverse();
						MATRIX __c_dU_inverse= - __c.inverse()*__c_dU*__c.inverse();
						double __j_dU=__j*(-__f.inverse().block(0,jj,DIM,1)*psi.transpose()).trace();
					
						SIGMAs_dU = mu_s*(2.0/3.0)*pow(__j,-1.0/3.0)*__j_dU*(MATRIX::Identity()-1.0/3.0*__c.trace()*__c.inverse());
						SIGMAs_dU +=mu_s*pow(__j,+2.0/3.0)*(-1.0/3.0*__c_dU.trace()*__c.inverse());
						SIGMAs_dU +=mu_s*pow(__j,+2.0/3.0)*(-1.0/3.0*__c.trace()*__c_dU_inverse);
						SIGMAs_dU +=0.5*kapa_s*(-2.0)*pow(__j,-3.0)*__j_dU*__c.inverse();
						SIGMAs_dU +=0.5*kapa_s*(pow(__j,-2.0)-1.0)*__c_dU_inverse;
				}
	  	else abort();
				sigmas_dU[jj] = __j*(-__f.inverse().block(0,jj,DIM,1)*psi.transpose()).trace() *__f.inverse()*SIGMAs*__f.inverse().transpose()
											+__j *(__f.inverse().block(0,jj,DIM,1)*psi.transpose()*__f.inverse())*SIGMAs*__f.inverse().transpose()
											+__j *__f.inverse() *SIGMAs*(__f.inverse().block(0,jj,DIM,1)*psi.transpose()*__f.inverse()).transpose()
											+__j *__f.inverse()*SIGMAs_dU*__f.inverse().transpose();
						
			}   
	 } 
  }
  




  template<int DIM>
  void SolidEuler<DIM>::MatrixBlock(EntryMatrix& A, const FemFunction& U, const FemFunction& N) const
  {;
    for (int j=0; j<N.size();++j)  // trial
      {
#define M N[j]
	point_M(j,U,M);
	
	for (int i=0; i<N.size();++i) // test
	  {
	    A.SetDofIndex(i,j);
	    Matrix(A,U,M , N[i]);
	  }
#undef M
      }
  }
  
  template class SolidEuler<3>;
  

}

/*-----------------------------------------*/
