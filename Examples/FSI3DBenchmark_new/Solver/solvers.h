/*----------------------------   solvers.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solvers_H
#define __solvers_H
/*----------------------------   solvers.h     ---------------------------*/

#include "sparse_umf.h"
#include "stdmultilevelsolver.h"
#include "stdsolver.h"
#include "stdsolver.h"
#include "alebasediscretization.h"
//#include "umfilu.h"
//#include "mumpsilu.h"
#include "fmatrixblock.h"
#include "lpsequation.h"
#include "facediscretization.h"


using namespace std;

namespace Gascoigne
{

  template<int DIM>
  class FSISolver : public StdSolver
  {
  private:
      
    DiscretizationInterface* NewDiscretization(int dimension, const std::string& discname);

     IluInterface *__IF,*__IS;

    //SparseUmf<FMatrixBlock<DIM+1> >         _LAP_M;    
    //SparseBlockMatrix<FMatrixBlock<DIM+1> > _LAP_A;
	
	std::string SolverLabel;
      
      //MumpsIlu<DIM+1>                         _LAP_M;

    public:
    
    DiscretizationInterface *&GetDiscretizationPointerPublic()
    {
      return _ZP;
    }
    FaceDiscretization *&GetFaceDiscretizationPointerPublic()
    {
      return _FZP;
    }
    
	  void SetSolverLabel(std::string solverlabel)
	  {
	  	SolverLabel=solverlabel;
	  }
	  std::string GetSolverLabel() const {return SolverLabel;}
	  
      vector<vector<double> >_precond;
      
      void do_precondition(VectorInterface gx) const;
      void undo_precondition(VectorInterface gx) const;
      
      FSISolver();

      virtual void ReInitMatrix() ;
  

     void SolveExtension(VectorInterface& x);


    void ComputeIlu(const VectorInterface& gu) const;
    void modify_ilu(IluInterface& I,int ncomp) const;

    void smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const ;
    void smooth_exact(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const;
      

      
    std::string GetName() const {return "FSI Solver";}
      
    AleBaseDiscretization* GetAleDiscretization() 
    {
      assert(dynamic_cast<AleBaseDiscretization*> (GetDiscretization()));
      return dynamic_cast<AleBaseDiscretization*> (GetDiscretization());
    }
    const AleBaseDiscretization* GetAleDiscretization()  const
    {
      assert(dynamic_cast<const AleBaseDiscretization*> (GetDiscretization()));
      return dynamic_cast<const AleBaseDiscretization*> (GetDiscretization());
    }
      
    void ReInitInterface(AleBaseDiscretization* ALEDISC);
    void reinit_element(int en, const nvector<int>& indices, 
			HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
			set<int>& fluid_nodes, set<int>& solid_nodes,int material);

    void reinit_interface_element(int en, const nvector<int>& indices, 
			  HASHMAP<int, std::vector<int> >& solid_interface_cells, 
			  HASHMAP<int, std::vector<int> >& fluid_interface_cells,
			  HASHSET<int> & interface_nodes,int material);

      

    void Form(VectorInterface& y, const VectorInterface& x, double d) const;
	
    void NewMesh(const MeshInterface* mp);
    void SetBoundaryVectorZero(VectorInterface& gf) const;
    void SetBoundaryVector(VectorInterface& gf) const;
    void DeleteSolidPressure(VectorInterface& gf) const;
      
    void AssembleMatrix(const VectorInterface& gu, double d);
      
    void PointVisu(const string& name, const GlobalVector& u, int iter) const;

        
    //      void ComputeIlu(const VectorInterface& gu) const;
      

      
    MatrixInterface* NewMatrix(int ncomp, const std::string& matrixtype);
    IluInterface* NewIlu(int ncomp, const std::string& matrixtype);
    
    void AddNodeVector(const std::string &name, const GlobalVector &q)
    {

      GetDiscretization()->AddNodeVector(name, &q);
    }  
  };
  
  

  
}

/*----------------------------   solvers.h     ---------------------------*/
/* end of #ifndef __solvers_H */
#endif
/*----------------------------   solvers.h     ---------------------------*/