#ifndef __solvers_h
#define __solvers_h

#include "stdmultilevelsolver.h"
#include "stdsolver.h"
#include "fmatrixblock.h"
//#include "sparse_umf.h"

namespace Gascoigne
{

  class MySolver : public StdSolver
  {
    //    SparseBlockMatrix<FMatrixBlock<2> >  DM;
    //    SparseUmf<FMatrixBlock<2> >          DUMF;
    
    GlobalVector                         DRHS;
    
  public:

    //    const SparseBlockMatrix<FMatrixBlock<2> >& GetDM() const   { return DM; }
    //    SparseBlockMatrix<FMatrixBlock<2> >& GetDM()               { return DM; }

    const GlobalVector& GetDRHS() const   { return DRHS; }
    GlobalVector& GetDRHS()    { return DRHS; }


    //    void SolveDual(VectorInterface& gz, VectorInterface& gf);
    //    void AddDualMatrix(const VectorInterface& gu, double d);
    //    void ReInitMatrix();

    void AssembleMatrix(const VectorInterface& gu, double d);
    void Form(VectorInterface& gy, const VectorInterface& gx, double d) const;      

  };
  

  class MyMLS : public StdMultiLevelSolver
  {
  public:
    
    SolverInterface* NewSolver(int solverlevel)
    {return new MySolver; }
  };


}

#endif
