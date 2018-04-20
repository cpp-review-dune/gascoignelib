#include "solvers.h"

#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "fmatrixblock.h"
#include <eigen3/Eigen/Dense>
#include "sparseblockmatrix.h"
#include "columndiagstencil.h"
#include <time.h> 

#include  "transformation3d.h"
#include  "baseq23d.h"
#include "gascoignehash.h"

using namespace std;



namespace Gascoigne
{


  template<int DIM>
  void FSISolver<DIM>::smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const 
  {			
    GlobalVector& gvX  = GetGV(x);
    const GlobalVector& gvY = GetGV(y);
    GlobalVector& gvH = GetGV(h);
    
    
    if (GetSolverData().GetLinearSmooth()=="vanka_jacobi"||GetSolverData().GetLinearSmooth()=="vanka_jacobi_av"||GetSolverData().GetLinearSmooth()=="vanka_gs")
      {

					string vanka_type=GetSolverData().GetLinearSmooth();
					double omega = GetSolverData().GetOmega();
				
					//const GascoigneMesh* M = dynamic_cast<const GascoigneMesh*> (GetMesh());
					//assert(M);
					
					const SparseBlockMatrix<FMatrixBlock<DIM> >* A = 
						dynamic_cast<const SparseBlockMatrix<FMatrixBlock<DIM> >*> (GetMatrix());
					assert(A);
					const ColumnDiagStencil* SA = dynamic_cast<const ColumnDiagStencil*> (A->GetStencil());
					assert(SA);
			    

					std::vector<int> smooth_weight;
					if(vanka_type=="vanka_jacobi" ||vanka_type=="vanka_jacobi_av")
					 	smooth_weight.resize(gvX.n(),0);
					 	
					GlobalVector smooth_zwisch;
					if(vanka_type=="vanka_jacobi_av")
						{
							smooth_zwisch.ncomp()=gvX.ncomp();
							smooth_zwisch.reservesize(gvX.n());
							smooth_zwisch.zero();
						}

			   	MatrixResidual(h,x,y);
					for(int iter=0; iter<niter; iter++)
						{
						  if((vanka_type=="vanka_jacobi" ||vanka_type=="vanka_jacobi_av")&&niter>0)
						  {
						  		MatrixResidual(h,x,y);
						  		for(vector<int>::iterator it_weight =smooth_weight.begin();it_weight !=smooth_weight.end();it_weight++)
						  			*it_weight=0;
						  		
						  }
						  if(vanka_type=="vanka_jacobi_av")
						  {  								
						  	smooth_zwisch.zero();	
						  }
							for (int p=0;p<Vanka_MV->patchindices_size();p++)
								{
						    
						    const nvector<int>& iop = Vanka_MV->Getpatchindices(p);
								//const HASHMAP<int,int>& inP=Vanka_MV->GetINP(p);
								//for (int i=0;i<iop.size();++i)
								//	inP[iop[i]]=i;

								int N = iop.size();
								//if(DIM==2) assert(N==9);
				        //if(DIM==3) assert(N==27);
								int ncomp = gvX.ncomp();
								// copy matrix & Vector

								Eigen::VectorXd H(N*ncomp);
					
								for (int r=0;r<N;++r)
									{
										int row = iop[r];
										// vector residuum
										for (int cr=0;cr<ncomp;++cr)
											H(ncomp*r+cr,0) = gvH(row,cr);			
									 }
								// local solve
//								Eigen::VectorXd X = Vanka_MV->GetMatrix(p).inverse()*H;
								Eigen::VectorXd X = Vanka_MV->GetLU(p).solve(H);

								// Update
								if(vanka_type=="vanka_gs")
								{
									for (int r=0;r<N;++r)
												for (int cr=0;cr<ncomp;++cr)
													gvX(iop[r],cr) += omega *  X(r*ncomp+cr,0);
		                
								}
								else if(vanka_type=="vanka_jacobi")
								{
									for (int r=0;r<N;++r)
									{
										if(smooth_weight[iop[r]]==0)
												for (int cr=0;cr<ncomp;++cr)
												{
														gvX(iop[r],cr) += omega *  X(r*ncomp+cr,0);
														smooth_weight[iop[r]]+=1;
												}
									}			
							  }
							  else if(vanka_type=="vanka_jacobi_av")
								{
									 for (int r=0;r<N;++r)
										{							
													for (int cr=0;cr<ncomp;++cr)
													{
														 smooth_zwisch(iop[r],cr) += omega *  X(r*ncomp+cr,0);
													}
											smooth_weight[iop[r]]++;		
										}		
								}
							
							
							
							
							//// Residual
								if(vanka_type=="vanka_gs") 
								{

										//MatrixResidual(h,x,y);
										A->MatrixResidualSome(Vanka_MV->Get_Vector_Columns(p), gvH,gvX,gvY);
    
/*											 for (auto row : Vanka_MV->Get_Vector_Columns(p))
											 {
												gvH.equ_node(row,row,gvY);												
												for (int pos=SA->start(row);pos<SA->stop(row);++pos)
												{				
													const int col = SA->col(pos); 

													const FMatrixBlock<DIM>& B = (*A->mat(pos));													

													for (int cr=0;cr<ncomp;++cr)
													  for (int cc=0;cc<ncomp;++cc)
													  	gvH(row,cr) -= B(cr,cc) * gvX(col,cc);
													  	
  										  	continue;
     												
											    
													if (inP.find(col)==inP.end()) continue;
													const int c=inP.find(col)->second;
													
													const FMatrixBlock<DIM>& B = (*A->mat(pos));
									
												  // hp[row] = hp[row] - omega * B * X[col]
													for (int cr=0;cr<ncomp;++cr)
													  for (int cc=0;cc<ncomp;++cc)
													  	gvH(row,cr) -= omega * B(cr,cc) * X(c*ncomp+cc);
													
													double* hp= &gvH(*row,0);

													double* Xp= &X(c*ncomp,0);
													for (int cr=0;cr<ncomp;++cr)
														{
															for (int cc=0;cc<ncomp;++cc)
																{
																			//*hp=GetGV(h)(*row,cr);
																			// *Xp=X(c*ncomp+cc,0)
																 			*hp-=omega *B(cr,cc)* *(Xp+cc);
																}
																++hp;
														 }

												} 
											 }
											 */
								}
					}
					if(vanka_type=="vanka_jacobi_av")
					{
						for (int i=0;i<gvX.n();++i)
							for (int ii=0;ii<gvX.ncomp();++ii)
								gvX(i,ii) +=smooth_zwisch(i,ii)/smooth_weight[i];
					}
			}

	  }
    else
      StdSolver::smooth(niter,x,y,h);
  
  
  }
  
  /* -------------------------------------------------------*/
  
  template<int DIM>
  void FSISolver<DIM>::invert_local_matrices() const
	{     

		const SparseBlockMatrix<FMatrixBlock<3> >* A = 
		dynamic_cast<const SparseBlockMatrix<FMatrixBlock<3> >*> (GetMatrix());
		assert(A);
	  
	  cout<<"Vanka_MV->ReInit start"<<endl;
    Vanka_MV->ReInit(A);
    cout<<"Vanka_MV->ReInit end"<<endl;
		/*			const GascoigneMesh* M = dynamic_cast<const GascoigneMesh*> (GetMesh());
					assert(M);
					const SparseBlockMatrix<FMatrixBlock<DIM> >* A = 
						dynamic_cast<const SparseBlockMatrix<FMatrixBlock<DIM> >*> (GetMatrix());
					assert(A);
					const ColumnDiagStencil* SA = dynamic_cast<const ColumnDiagStencil*> (A->GetStencil());
					assert(SA);
				
				 Vanka_MV->delete_Vector_Columns();
							
				 for (int p=0;p<NumberofSmoothingCells();++p)
					{
						  const nvector<int>& iop = IndicesSmoothingCell(p);
						
					std::map<int,int> inP;
					for (int i=0;i<iop.size();++i)
						inP[iop[i]]=i;

					int N = iop.size();
					//if(DIM==2) assert(N==9);
				  //if(DIM==3) assert(N==27);
	 				
					int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();

					// copy matrix & Vector
					Eigen::MatrixXd P(N*ncomp,N*ncomp);
					

					for (int r=0;r<N;++r)
						{					
						 
							int row = iop[r];
							// matrix
							for (int pos=SA->start(row);pos<SA->stop(row);++pos)
								{
									int col = SA->col(pos);
									if(GetSolverData().GetLinearSmooth()=="vanka_gs") 
										Vanka_MV->insert_Vector_Columns(p,col);
									
									if (inP.find(col)==inP.end()) continue;
									int c = inP[col];
									const FMatrixBlock<DIM>& B = (*A->mat(pos));
									for (int cr=0;cr<ncomp;++cr)
										for (int cc=0;cc<ncomp;++cc)
											P(ncomp*r+cr,ncomp*c+cc) = B(cr,cc);
								}

						 }
					// local solve
				  Vanka_MV->insert(p,P.inverse());
				  
				  
			}*/
	}
/* -------------------------------------------------------*/
	
  template<int DIM>
  void FSISolver<DIM>::ComputeIlu() const
{
		if ( (GetSolverData().GetLinearSmooth()=="vanka_jacobi"||GetSolverData().GetLinearSmooth()=="vanka_jacobi_av"||GetSolverData().GetLinearSmooth()=="vanka_gs" )&& _directsolver==false)
    { 
    	invert_local_matrices();
    	cout<<"??????????????????????????"<<endl;
    }
    else 
    StdSolver::ComputeIlu() ;
}

/* -------------------------------------------------------*/

  template<int DIM>
  void FSISolver<DIM>::ComputeIlu(const VectorInterface& gu) const
{
		if ( (GetSolverData().GetLinearSmooth()=="vanka_jacobi"||GetSolverData().GetLinearSmooth()=="vanka_jacobi_av"||GetSolverData().GetLinearSmooth()=="vanka_gs" )&& _directsolver==false)
    	{
    		invert_local_matrices();
    	}
    else
    	StdSolver::ComputeIlu(gu);
}
/* -------------------------------------------------------*/
  template<int DIM>
  void FSISolver<DIM>::RegisterMatrix()
{
	StdSolver::RegisterMatrix();

		if ( (GetSolverData().GetLinearSmooth()=="vanka_jacobi"||GetSolverData().GetLinearSmooth()=="vanka_jacobi_av"||GetSolverData().GetLinearSmooth()=="vanka_gs" )&& _directsolver==false)
  {
  	const GascoigneMesh* M = dynamic_cast<const GascoigneMesh*> (GetMesh());
		assert(M);
		cout<<"RegisterMatrix"<<endl;
		if(Vanka_MV==NULL)
			Vanka_MV=new Vanka_Matrix_Vector;
			
		//Vanka_MV->resize_Matrix_Vector(NumberofSmoothingCells());
		Vanka_MV->build_patches(M);
  }	
}
/* -------------------------------------------------------*/
  template<int DIM>
   FSISolver<DIM>::~FSISolver()
	{
		if ( (GetSolverData().GetLinearSmooth()=="vanka_jacobi"||GetSolverData().GetLinearSmooth()=="vanka_jacobi_av"||GetSolverData().GetLinearSmooth()=="vanka_gs" )&& _directsolver==false)	
				Vanka_MV->clear();
				
	}
	
	
  template class FSISolver<2>;
  template class FSISolver<3>;
  template class FSIMultiLevelSolver<2>;
  template class FSIMultiLevelSolver<3>;
}
