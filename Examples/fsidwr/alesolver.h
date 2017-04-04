/*----------------------------   alesolver.h     ---------------------------*/
/*      $Id: alesolver.h,v 1.11 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __alesolver_H
#define __alesolver_H
/*----------------------------   alesolver.h     ---------------------------*/


#include "stdsolver.h"
#include "alebasediscretization.h"
#include "pointmatrix.h"
#include "sparseblockmatrix.h"

namespace Gascoigne
{

  class AleSolver : public StdSolver
  {
  private:

    DiscretizationInterface* NewDiscretization(int dimension, const std::string& discname);

    

    // Submatrices for Fluid & Solid
    // used for smoothing in a domain-partitioning approach
    bool __DUALMATRIX;
    bool __splittingsmoother;
    bool __splitting_fluid_exact;
    bool __splitting_solid_exact;
    
    bool cond_computed;
    MatrixInterface *__AF,*__AS;
    IluInterface    *__IF,*__IS;

    // local vectors in fluid- and solid-domain for linear solver
    mutable GlobalVector __Xf,__Xs, __Hf, __Hs;
    
    
  public:
    AleSolver();
    ~AleSolver();
    int  Richardson(const MatrixInterface &A, 
		    GlobalVector &x, const GlobalVector &b,
		    const IluInterface &M, 
		    int &max_iter, 
		    double &tol) const ;
    

    int  BiCGSTAB(const MatrixInterface &A, 
		  GlobalVector &x, const GlobalVector &b,
		  const IluInterface &M, 
		  int &max_iter, 
		  double &tol) const;

    bool SplittingSmoother() const   { return __splittingsmoother; }
    bool SplittingSolidExact() const { return __splitting_solid_exact; }
    bool SplittingFluidExact() const { return __splitting_fluid_exact; }



    ////////// SPLITTING SOLVER
    void Fluidg2l(GlobalVector& Xf, const GlobalVector& x) const;
    void Solidg2l(GlobalVector& Xs, const GlobalVector& x) const;
    void Fluidl2g(GlobalVector& x, double s, const GlobalVector& Xf) const;
    void Solidl2g(GlobalVector& x, double s, const GlobalVector& Xs) const;
    void FluidInterfaceZero(GlobalVector& X, int c) const;
    void SolidInterfaceZero(GlobalVector& X, int c) const;
    
    void SolveFluidExact(VectorInterface& x,const VectorInterface& h) const;
    void SolveSolidExact(VectorInterface& x,const VectorInterface& h) const;
    
    void SolveFluidIterative(int niter, double TOL, VectorInterface& x,const VectorInterface& h) const;
    void SolveSolidIterative(int niter, double TOL, VectorInterface& x,const VectorInterface& h) const;


    
    void Visu(const std::string& name, const VectorInterface& gu, int i) const;
    
    std::string GetName() const {return "ALE Solver";}

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
    
    const HASHSET<int>& GetInterfaceNodes() const { return GetAleDiscretization()->GetInterfaceNodes(); }
    const std::vector<int>&   GetFluidL2G() const
    { return GetAleDiscretization()->GetFluidL2G(); }
    const std::vector<int>&   GetSolidL2G() const
    { return GetAleDiscretization()->GetSolidL2G(); }


    void ReInitInterface(AleBaseDiscretization* ALEDISC);
    void reinit_element_2d(int en, const nvector<int>& indices, 
			   HASHMAP<int, std::vector<int> >& solid_interface_cells, 
			   HASHMAP<int, std::vector<int> >& fluid_interface_cells,
			   HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
			   HASHSET<int> & interface_nodes,
			   std::set<int>& fluid_nodes, std::set<int>& solid_nodes);
    void reinit_element_3d(int en, const nvector<int>& indices, 
			   HASHMAP<int, std::vector<int> >& solid_interface_cells, 
			   HASHMAP<int, std::vector<int> >& fluid_interface_cells,
			   HASHSET<int> & fluid_cells, HASHSET<int> & solid_cells,
			   HASHSET<int> & interface_nodes,
			   std::set<int>& fluid_nodes, std::set<int>& solid_nodes);

    void NewMesh(const MeshInterface* mp);
    void SetBoundaryVectorZero(VectorInterface& gf) const;
    void SetBoundaryVector(VectorInterface& gf) const;
    
    void Form(VectorInterface& gy, const VectorInterface& gx, double d) const;
    void AssembleMatrix(const VectorInterface& gu, double d);
    void AssembleDualMatrix(const VectorInterface& gu, double d);
    void ComputeIlu(const VectorInterface& gu) const;
    MatrixInterface* NewMatrix(int ncomp, const std::string& matrixtype); 
    IluInterface* NewIlu(int ncomp, const std::string& matrixtype); 


    void smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const;

    template<class FF, class FS, class FA> void SplitMatrix(SparseBlockMatrix<FF>* AF,SparseBlockMatrix<FS>* AS, const SparseBlockMatrix<FA>* A);


    double LargestEV(MatrixInterface& M);
    double SmallestEV(MatrixInterface& M);
    void SwapUV(MatrixInterface& M);


    void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
    
  };

}



/*----------------------------   alesolver.h     ---------------------------*/
/* end of #ifndef __alesolver_H */
#endif
/*----------------------------   alesolver.h     ---------------------------*/
