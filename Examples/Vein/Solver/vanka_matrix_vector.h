/*----------------------------  vanka_matrix_vector.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __vanka_matrix_vector_H
#define __vanka_matrix_vector_H
/*----------------------------   vanka_matrix_vector.h     ---------------------------*/



#include <eigen3/Eigen/Dense>
#include <iostream>
#include <map>
#include <set>
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "fmatrixblock.h"
#include <eigen3/Eigen/Dense>
#include "sparseblockmatrix.h"
#include "columndiagstencil.h"
#include "gascoignehash.h"
#include <fstream> 

#include "vanka_matrix_vector_base.h"


namespace Gascoigne
{

    class VertexCompare
    {
    public:
      bool operator()(const Vertex3d& v1, const Vertex3d& v2) const
    	{
    	  /*if (fabs(v1.y()-v2.y())<1.e-10)
    	  {
    	    return (v1.x()+1.e-10<v2.x());
    	  }
    	  return v1.y()+1.e-10<v2.y();
    	  */
    	  if (fabs(v1.z()-v2.z())<1.e-10)
    	  {
    	  	if (fabs(v1.y()-v2.y())<1.e-10)
    	  	   return (v1.x()+1.e-10<v2.x()); 
    	  	else 
    	    	return (v1.y()+1.e-10<v2.y());
    	  }
    	  return v1.z()+1.e-10<v2.z();
    	}
    };
    
    class VertexCompare_y
    {
    public:
      bool operator()(const Vertex3d& v1, const Vertex3d& v2) const
    	{
    	  	if (fabs(v1.y()-v2.y())<1.e-10)
    	  	   return (v1.x()+1.e-10<v2.x()); 
    	  	else 
    	    	return (v1.y()+1.e-10<v2.y());
    	}
    };
    
    class VertexCompare_x
    {
    public:
      bool operator()(const Vertex3d& v1, const Vertex3d& v2) const
    	{
    	//Es kann bei passieren das 
    		if(fabs(v1.x()-v2.x())<1.e-10)
    			return (v1.y()+1.e-10<v2.y());
    	  	else return (v1.x()+1.e-10<v2.x());  	
    	  	
    	}
    };
    
    template<int NCOMP>
    class Vanka_Matrix_Vector : public virtual Vanka_Matrix_Vector_base
    {
    private:
			//std::vector<Eigen::MatrixXd>	Inverse_Matrix_Vector;
   		std::vector<Eigen::FullPivLU<Eigen::MatrixXd> >	Inverse_LU_Vector;		
			std::vector<std::vector<int>> Vector_Columns_Patch;
			
	    std::vector<std::vector<int> > patchindices;
			std::vector<HASHMAP<int,int> > INP;	
    public:
    
     
    	void build_patches(const GascoigneMesh* M)
    	{
        	assert(M);
        	const GascoigneMesh3d* M3 = dynamic_cast<const GascoigneMesh3d*>(M);
        	assert(M3);
        	
          patchindices.clear();
        	
        	std::map<int,std::map<Vertex3d,int, VertexCompare>> X; //VertexHash, VertexEqual> X;
        	std::map<int,std::vector<int>> Z;
        	
        	//Zaehlen wieviele Patches in Urvaterblock
        	std::map<int,int> patchesofsamecol;
        	for (int p=0;p<M3->npatches();++p)
        		{
        		 if(patchesofsamecol.find(M3->material_Vanka_patch(p))!=patchesofsamecol.end())
        		 	patchesofsamecol[M3->material_Vanka_patch(p)]++;
        		 else	
        			patchesofsamecol[M3->material_Vanka_patch(p)]=1;
        		}
        	
        	//Schleife uber alle Patches und Blocken der Patche mit selben Urvater
        	//Sortieren der Patche jeweils entsprechend dem lokalen Koordinatensystem des Urvaters	
        	for (int p=0;p<M3->npatches();++p)
        	{
        	  const IntVector& iop = *(M3->IndicesOfPatch(p));
        	  //Umrechnung der Koordinaten des Patchmittelpunktes in lokales Koordinatensystem
        	  Eigen::MatrixXd BasisM(3,3);
        	  Eigen::Vector3d NewBasis;
        	  for(int j=0;j<3;j++) 
        	  	NewBasis(j)=M3->vertex3d(iop[13])[j];
        	  for(int ii=0;ii<3;ii++)
        	  {
        	   for(int iii=0;iii<3;iii++)
        	  		{ 
        	  			BasisM(iii,ii)=M3->Vanka_basis_patch(p)[ii][iii];
        	  		}
        	  }
						
        	  NewBasis=	BasisM.inverse()*NewBasis;
		
						//Vertex3d mid(M3->vertex3d(iop[13]).x(),M3->vertex3d(iop[13]).y(),M3->vertex3d(iop[13]).z());
						Vertex3d mid(NewBasis(0),NewBasis(1),NewBasis(2));
						
        		X[M3->material_Vanka_patch(p)][mid]=p;
        		//patchindices.push_back(iop);
        		
        	}
        	
        	//Problem bei komplexeren Geometrien Patche in X nur in z-Koordinate sortiert
        	//Blocken von jeweils n*n patches und sortieren dieser Bloecke in y Richtung
        	//Blocken von n patches in jedem Block und sortieren in x-Richtung
        	for (auto iterX : X)
        	{
        	  //Schleife uber alle Urvaterbloecke
        	  //iterX.second Mittelpunkte sortiert nach z

						int n=round(pow(patchesofsamecol[iterX.first],1.0/3.0));

        	  //Fasse n*n Patches(sortiert nach z) zusammen und sortiere diese nach y
        	  std::vector<std::map<Vertex3d,int,VertexCompare_y>> Blocks;
        	  Blocks.resize(n);
        	  int block_counter=0;

        	  for (auto iter_iter_X : iterX.second)
        	  {
        	   //Schleife uber alle Patches in Urvaterzelle
        	   Vertex3d mid=iter_iter_X.first;
        	   int p=iter_iter_X.second;
        	  
        	   int block=block_counter/(n*n);
        	   Blocks[block][mid]=p;
        	   block_counter++;
        		}
        		
						//Fasse n Patches(sortiert nach y) in jedem der n Bloecke zusammen und sortiere diese nach x
        	  for (auto iter_Blocks : Blocks)
        	  {
        	    //iter_Blocks Mittelpunkte eines Blocks sortiert nach y
							std::vector<std::map<Vertex3d,int,VertexCompare_x>> Sub_Blocks;
							Sub_Blocks.resize(n);
							
		      	  int sub_block_counter=0;
		      	  
		      	  for(auto iter_patches_in_Block : iter_Blocks)
				    	 {
				    	   //Schleife uber alle Patches in Block
						  	 Vertex3d mid=iter_patches_in_Block.first;
						  	 int p=iter_patches_in_Block.second;
						  	 int block2=sub_block_counter/n;
						  	 Sub_Blocks[block2][mid]=p;

						  	 sub_block_counter++;
				    		}		
						//Speichere die Patches entsprechend ihrer Sortierung ab	
							for(auto iter_Sub_Blocks : Sub_Blocks)
							 for(auto iter_patch: iter_Sub_Blocks) 
							 	Z[iterX.first].push_back(iter_patch.second);
							 	
							 	
        		}

        	//	patchindices.push_back(iop);
        	}
        	//int combine_n_patches=1;
        	//cout<<"patchesofsamecol[0]"<<patchesofsamecol[0]<<endl;

        	
					int count_iterbox=0;
        	
        	for(auto iter_Box : Z)
        	{
        	
        	  int combine_n_patches;
        	  if(patchesofsamecol[iter_Box.first]==1) combine_n_patches=1;
        	  if(patchesofsamecol[iter_Box.first]==8) combine_n_patches=2;
        	  if(patchesofsamecol[iter_Box.first]==8*8) combine_n_patches=4;
        	  if(patchesofsamecol[iter_Box.first]==8*8*8) combine_n_patches=8;
        	  
        	  int n_smooth_patche=iter_Box.second.size()/combine_n_patches;

		      	/*std::vector<int> X_Vec;
		      	for (auto iterX : iter_Box.second)
		      	{
		      		X_Vec.push_back(iterX.second);
		      	}	
		      	*/
		        for (int i=0;i<n_smooth_patche;i++)
		        {

		          std::set<int> tmp_set;
		          for(int ii=0;ii<combine_n_patches;ii++)
		          {
		             const IntVector& iop = *(M3->IndicesOfPatch(iter_Box.second[i*combine_n_patches+ii]));
		           	 tmp_set.insert(iop.begin(),iop.end()); 		
		          }
		          std::vector<int> tmp_vec(tmp_set.begin(),tmp_set.end());
		          //for(auto iii : tmp_set)
		          //  tmp_vec.push_back(iii);
		           
		          patchindices.push_back(tmp_vec);
		          	
		        }  	
          }  
           /* vector<int> tmppatch;
            for (auto ij : tmp)
            	for (auto ik : ij.second)
                tmppatch.push_back(ik);
                        
            patchindices.push_back(tmppatch);
//            cout << tmppatch.size() << endl;
          }*/
              
          INP.resize(patchindices.size());
          for (int p=0;p<patchindices.size();++p)
          {
          	INP[p].clear();
          	for (int i=0;i<patchindices[p].size();++i)
          		INP[p][patchindices[p][i]]=i;
          }
          
         /* ofstream func_log("blub.txt");
    			func_log.precision(12);
          for (int p=0;p<patchindices.size();++p)
           {
           	for(auto patches : patchindices[p])
           	{
           		func_log <<M3->vertex3d(patches)[0]<<" "<<M3->vertex3d(patches)[1]<<" "<<M3->vertex3d(patches)[2]<<" ";
           	}
           	func_log <<endl;
           }
           */
    	}
    

       void ReInit(const MatrixInterface * A_Matrix) 
      {
      	
      	const SparseBlockMatrix<FMatrixBlock<NCOMP> >* A = 
		dynamic_cast<const SparseBlockMatrix<FMatrixBlock<NCOMP> >*> (A_Matrix);
		assert(A);
		
        assert(A);
      	const ColumnDiagStencil* SA = dynamic_cast<const ColumnDiagStencil*> (A->GetStencil());
    		assert(SA);

			  int npatch = patchindices.size();
				//Inverse_Matrix_Vector.resize(npatch);
				Inverse_LU_Vector.clear();
   			Vector_Columns_Patch.resize(npatch);
				Eigen::MatrixXd Matrix_on_Block;
				
				for (int p=0;p<npatch;++p)
				{
					int np = patchindices[p].size();
					
					Matrix_on_Block.resize(np*NCOMP,np*NCOMP);
					Matrix_on_Block.setZero();
					

					Vector_Columns_Patch[p].clear();
					set<int> columnsinpatch;	
					
					for (int r=0;r<np;++r)
						{
							int row = patchindices[p][r];
							// matrix
							for (int pos=SA->start(row);pos<SA->stop(row);++pos)
								{
									int col = SA->col(pos);
									columnsinpatch.insert(col);
									
									if (INP[p].find(col)==INP[p].end()) continue;
									int c = INP[p][col];
									const FMatrixBlock<NCOMP>& B = (*A->mat(pos));
									for (int cr=0;cr<NCOMP;++cr)
										for (int cc=0;cc<NCOMP;++cc)
											Matrix_on_Block(NCOMP*r+cr,NCOMP*c+cc) = B(cr,cc);
								}
						 }
				  Inverse_LU_Vector.push_back(Eigen::FullPivLU<Eigen::MatrixXd> (Matrix_on_Block ) );

//  				 Inverse_Matrix_Vector[p]=Inverse_Matrix_Vector[p].inverse();
//             Inverse_LU_Vector[p]=Inverse_Matrix_Vector[p].inverse();

 		       for (auto it : columnsinpatch)
		       	 Vector_Columns_Patch[p].push_back(it);
				}
				
      }
      //Eigen::MatrixXd& GetMatrix(int n) 
      //{
      //	return Inverse_Matrix_Vector[n];
      //}
      Eigen::FullPivLU<Eigen::MatrixXd>& GetLU(int n) 
      {
      	return Inverse_LU_Vector[n];
      }
      //void delete_Matrix_Vector()  
      //{
      // Inverse_Matrix_Vector.clear();
      //}
     //void resize_Matrix_Vector(int n)  
      //{
      // Inverse_Matrix_Vector.resize(n);
      //}
      
      int patchindices_size()
      {
       return patchindices.size();
      }
      
      const std::vector<int>&  Getpatchindices( int p)
      {
       return patchindices[p];
      }
      
      const HASHMAP<int,int>&  GetINP( int p)
      {
       return INP[p];
      }
      
      const std::vector<int>& Get_Vector_Columns(int p)
      {
        return Vector_Columns_Patch[p];
      }
      
      std::vector<int>::iterator Get_Vector_Columns_begin(int n)
      {
        return Vector_Columns_Patch[n].begin();
      }
      
      std::vector<int>::iterator Get_Vector_Columns_end(int n)
      {
        return Vector_Columns_Patch[n].end();
      }
      void clear()  
      {
      	Vector_Columns_Patch.clear();
      	Inverse_LU_Vector.clear();		
	      patchindices.clear();
			  INP.clear();
      }
  
    };
  }
  


/*----------------------------   vanka_matrix_vector.h     ---------------------------*/
/* end of #ifndef __vanka_matrix_vector_H */
#endif
/*----------------------------   vanka_matrix_vector.h     ---------------------------*/
