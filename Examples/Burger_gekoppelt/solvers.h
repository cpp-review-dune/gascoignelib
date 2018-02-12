#ifndef __solvers_h
#define __solvers_h

#include "stdmultilevelsolver.h"
#include "stdsolver.h"
#include "fmatrixblock.h"
#include "sparse_umf.h"

namespace Gascoigne
{

  class MySolver : public StdSolver
  {
    SparseBlockMatrix<FMatrixBlock<2> >  DM;
    SparseUmf<FMatrixBlock<2> >  DUMF;
    
    GlobalVector                         DRHS;
    
  public:

    const SparseBlockMatrix<FMatrixBlock<2> >& GetDM() const   { return DM; }
    SparseBlockMatrix<FMatrixBlock<2> >& GetDM()               { return DM; }

    const GlobalVector& GetDRHS() const   { return DRHS; }
    GlobalVector& GetDRHS()    { return DRHS; }


    void SolveDual(VectorInterface& gz, VectorInterface& gf)
    {
      GlobalVector& z = GetGV(gz);

      DM.transpose();
      z.zero();

      GetGV(gf) = DRHS;
      SetBoundaryVectorZero(gf);

      vector<int> perm(z.n());
      for (int i=0;i<z.n();++i) perm[i]=i;
      DUMF.ConstructStructure(perm, DM);
      DUMF.copy_entries(&DM);
      DUMF.compute_ilu();

      z=DRHS;
      DUMF.solve(z);
      
    }

    
    void AddDualMatrix(const VectorInterface& gu, double d)
    {
      const GlobalVector& u = GetGV(gu);
      
      HNAverage(gu);
      HNAverageData();

      GetDiscretization()->Matrix(DM,u,*GetProblemDescriptor()->GetEquation(),d);
            
      DirichletMatrix();
      HNZero(gu);
      HNZeroData();

      ////
      const DomainRightHandSide *RHS = dynamic_cast<const DomainRightHandSide *>(GetProblemDescriptor()->GetRightHandSide());
      assert(RHS);
      GetDiscretization()->AddNodeVector("U",&u);
      GetDiscretization()->Rhs(DRHS,*RHS, d);
      GetDiscretization()->DeleteNodeVector("U");
	    
    }
    
    void ReInitMatrix()  
    {
      // hier aus stdsolver
      GetDiscretization()->InitFilter(GetPfilter());
      SparseStructure SA;
      GetDiscretization()->Structure(&SA);
      
      if (GetFaceDiscretization())
	GetFaceDiscretization()->Structure(&SA);
      
      AddPeriodicNodes(&SA);
      
      GetMatrix()->ReInit(&SA);
      GetIlu()->ReInit(&SA);

      // hier neu
      DM.ReInit(&SA);
      DUMF.SetMatrix(&DM);
      DUMF.ReInit(&SA);
    }

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
