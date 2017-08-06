/*----------------------------   solvers.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solvers_H
#define __solvers_H
/*----------------------------   solvers.h     ---------------------------*/


#include "stdmultilevelsolver.h"
#include "stdsolver.h"
#include "alebasediscretization.h"
#include "umfilu.h"

using namespace std;

namespace Gascoigne
{

  template<int DIM>
    class FSISolver : public StdSolver
    {
    private:
      
      DiscretizationInterface* NewDiscretization(int dimension, const std::string& discname);

      IluInterface *__IF,*__IS;
      
      
    public:


      vector<vector<double> >_precond;
      
      void do_precondition(VectorInterface gx) const;
      void undo_precondition(VectorInterface gx) const;
      
      FSISolver();

      // >>>>>>>>>>>>>>>>> ILU STUFF
      void ReInitMatrix() ;
      void ComputeIlu(const VectorInterface& gu) const;
      void smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const ;
      void smooth_exact(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const;
      
      void  Precondition(const IluInterface& M, GlobalVector &x) const;
      int  BiCGSTAB(const MatrixInterface &A, 
		    GlobalVector &x, const GlobalVector &b,
		    const IluInterface &M, 
		    int &max_iter, double &tol) const ;
      void  Precondition_SINGLE(const IluInterface& M, GlobalVector &x) const;
      int  BiCGSTAB_SINGLE(const MatrixInterface &A, 
			   GlobalVector &x, const GlobalVector &b,
			   const IluInterface &M, 
			   int &max_iter, double &tol) const ;

    


      
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
			  HASHMAP<int, std::vector<int> >& solid_interface_cells, 
			  HASHMAP<int, std::vector<int> >& fluid_interface_cells,
			  HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
			  HASHSET<int> & interface_nodes,
			  set<int>& fluid_nodes, set<int>& solid_nodes);
      

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
      void ComputeIlu(VectorInterface& u);
    };
  
}

/*----------------------------   solvers.h     ---------------------------*/
/* end of #ifndef __solvers_H */
#endif
/*----------------------------   solvers.h     ---------------------------*/
