#include "solvers.h"

extern double CHIGAUSS;

namespace Gascoigne{

  /*
  void MySolver::SolveDual(VectorInterface& gz, VectorInterface& gf)
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
  
    
  void MySolver::AddDualMatrix(const VectorInterface& gu, double d)
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
  
  void MySolver::ReInitMatrix()  
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
  */
  
  void MySolver::AssembleMatrix(const VectorInterface& gu, double d)
  {
    CHIGAUSS = 0.5-0.5*sqrt(1.0/3.0);
    StdSolver::AssembleMatrix(gu,d*0.5);
    CHIGAUSS = 0.5+0.5*sqrt(1.0/3.0);
    StdSolver::AssembleMatrix(gu,d*0.5);
  }
  
void MySolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  CHIGAUSS = 0.5-0.5*sqrt(1.0/3.0);
  StdSolver::Form(gy,gx,d*0.5);
  CHIGAUSS = 0.5+0.5*sqrt(1.0/3.0);
  StdSolver::Form(gy,gx,d*0.5);
}
  
}
