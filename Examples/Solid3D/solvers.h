/*----------------------------   solvers.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solvers_H
#define __solvers_H
/*----------------------------   solvers.h     ---------------------------*/


#include "stdmultilevelsolver.h"
#include "stdsolver.h"
#include <eigen3/Eigen/Dense>
#include "vanka_matrix_vector.h"


using namespace std;

namespace Gascoigne
{

  template<int DIM>
    class FSISolver : public StdSolver
    {
    private:
      Vanka_Matrix_Vector* Vanka_MV=NULL;
    public:
			~FSISolver();
			
      void smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const ;
      void invert_local_matrices() const  ;
      void ComputeIlu(const VectorInterface& gu) const;
      void ComputeIlu() const;
      void RegisterMatrix();
      void Anisotropy() const ;
      void Transformation(FemInterface::Matrix& T, int iq) const;
      
      int NumberofSmoothingCells()const
      { 
      	
      	const GascoigneMesh* M = dynamic_cast<const GascoigneMesh*> (GetMesh());
      	return M->npatches();
      	//return 2*M->npatches();
				//return 8*M->nq4patches();
      }
      const nvector<int> IndicesSmoothingCell(int p)const
      {
        const GascoigneMesh* M = dynamic_cast<const GascoigneMesh*> (GetMesh());
        return *(M->IndicesOfPatch(p));
        /*    
            const GascoigneMesh* M = dynamic_cast<const GascoigneMesh*> (GetMesh());
					  nvector<int> iop;
						nvector<int> cellsofpatch = M->GetPatchIndexHandler().GetPatch2Cell(p/2);
						
						set<int> indices_half_patch_1;
						set<int> indices_half_patch_2;
						
						int pp1[]={0,1,4,5};
						set<int> patch_part_1(pp1,pp1+4);
						int pp2[]={2,3,6,7};
						set<int> patch_part_2(pp2,pp2+4);

						for(set<int>::iterator itttt=patch_part_1.begin();itttt!=patch_part_1.end();itttt++)
						{
							nvector<int> indcell=M->IndicesOfCell(cellsofpatch[*itttt]);
							vector<int>::iterator it;
							for (it=indcell.begin(); it!=indcell.end(); ++it)			
								indices_half_patch_1.insert(*it);
						}
						for(set<int>::iterator itttt=patch_part_2.begin();itttt!=patch_part_2.end();itttt++)
						{
							nvector<int> indcell=M->IndicesOfCell(cellsofpatch[*itttt]);
							vector<int>::iterator it;
							for (it=indcell.begin(); it!=indcell.end(); ++it)			
								indices_half_patch_2.insert(*it);
						}
					  
					if(p%2==0)
					 {
					 	 set<int>::iterator it;
					 	 for (it=indices_half_patch_1.begin(); it!=indices_half_patch_1.end(); ++it)			
							iop.push_back(*it);
					 }
					else	
					 {
					 	 set<int>::iterator it;
					 	 for (it=indices_half_patch_2.begin(); it!=indices_half_patch_2.end(); ++it)			
							iop.push_back(*it);
					 }
					 return iop;
					*/ 
	
/*
					 const GascoigneMesh* M = dynamic_cast<const GascoigneMesh*> (GetMesh());
					  nvector<int> iop;
						nvector<int> cellsofpatch = M->GetPatchIndexHandler().GetQ4Patch2Cell(p/8);
						set<int> indices_half_patch;
						
						vector<int> pp; pp.push_back(0);pp.push_back(1);pp.push_back(2);pp.push_back(3);pp.push_back(4);pp.push_back(5);pp.push_back(6);pp.push_back(7);

						
						for(int ii=0;ii<8;ii++)
						{
							nvector<int> indcell=M->IndicesOfCell(cellsofpatch[pp[p%8]*8+ii]);
							vector<int>::iterator it;
							for (it=indcell.begin(); it!=indcell.end(); ++it)			
								indices_half_patch.insert(*it);
						}
						
					 	 set<int>::iterator it;
					 	 for (it=indices_half_patch.begin(); it!=indices_half_patch.end(); ++it)			
							iop.push_back(*it);
					return iop;
			*/		 
    }  
    };
  
  

  template<int DIM>
    class FSIMultiLevelSolver : public StdMultiLevelSolver
    {
    public:
      
            
      std::string GetName() const {return "FSI MultiLevelSolver";}
	  
      SolverInterface* NewSolver(int solverlevel)
      { return new FSISolver<DIM>; }
      
      const FSISolver<DIM>* GetFSISolver(int l) const
      {
	assert(dynamic_cast<const FSISolver<DIM>* > (GetSolver(l)));
	return dynamic_cast<const FSISolver<DIM>* > (GetSolver(l));
      }
    };
  
}

/*----------------------------   solvers.h     ---------------------------*/
/* end of #ifndef __solvers_H */
#endif
/*----------------------------   solvers.h     ---------------------------*/
