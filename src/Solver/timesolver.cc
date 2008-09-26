#include "timesolver.h"
#include  "cg.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

void TimeSolver::SetTimeData(double _dt, double _theta, double _time) 
{
  dt    = _dt; 
  theta = _theta;
  time  = _time;
  assert(dt>0.);
  assert(theta>0.);
  GetProblemDescriptor()->SetTime(time,dt);
}

/*-----------------------------------------*/

void TimeSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->GetNcomp();
  
  if (GetMassMatrixPointer()==NULL)
    GetMassMatrixPointer() = NewMassMatrix(ncomp,_matrixtype);
    
  StdSolver::RegisterMatrix();
}

/*-----------------------------------------*/

void TimeSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  const Equation* EQ = PDX.GetEquation();
  if (EQ) 
    {
      _TP.reservesize(EQ->GetNcomp(),EQ->GetNcomp(),0.);
      EQ->SetTimePattern(_TP);
    }   
  StdSolver::SetProblem(PDX);
}

/*-----------------------------------------*/

void TimeSolver::ReInitMatrix() 
{
  GetDiscretization()->InitFilter(_PF);
  SparseStructure SA;
  GetDiscretization()->Structure(&SA);
  
  GetMatrix()->ReInit(&SA);
  GetIlu()   ->ReInit(&SA);
  
  GetMassMatrix()->ReInit(&SA);
  GetMassMatrix()->zero();
  GetDiscretization()->MassMatrix(*GetMassMatrix()); 
}

/*-----------------------------------------*/

void TimeSolver::AssembleMatrix(const VectorInterface& gu, double d)
{
  StdSolver::AssembleMatrix(gu,d);
  
  double scale = d/(dt*theta);
  GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),_TP,scale);

  StdSolver::DirichletMatrix();
}

/*-----------------------------------------*/
  
void TimeSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  StdSolver::Form(gy,gx,d);
  
  double scale = d/(dt*theta);
  MassMatrixVector(gy,gx,scale);
}

/*-----------------------------------------*/

void TimeSolver::MassMatrixVector(VectorInterface& gf, const VectorInterface& gu, double d) const
{
        GlobalVector& f = GetGV(gf);
  const GlobalVector& u = GetGV(gu);
  GetMassMatrix()->vmult_time(f,u,_TP,d);
}

/*-----------------------------------------*/

void TimeSolver::InverseMassMatrix(VectorInterface& u, const VectorInterface& f, CGInfo& info)
{
  CG<TimeSolver,VectorInterface> cg(*this);
  cg.solve(u,f,info);
}

/*-----------------------------------------*/

void TimeSolver::precondition(VectorInterface& u, const VectorInterface& f)
{
  Equ(u,1.,f);
  GetMassMatrix()->PrepareJacobi(1.);
  GetMassMatrix()->Jacobi(GetGV(u));
}

/*-----------------------------------------*/

void TimeSolver::cgvmult(VectorInterface& y, const VectorInterface& x, double d) const
{
  MassMatrixVector(y,x,d);
}

/*-----------------------------------------*/

}
