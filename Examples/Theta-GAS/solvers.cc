#include "solvers.h"
#include "problem.h"
#include "diracrighthandside.h"

namespace Gascoigne
{
  void ThetaMultiLevelSolver::SolveAdjoint(VectorInterface& z,  VectorInterface& f, CGInfo& info)
  {
    Zero(z);
    ////////////////// im dualen nur dirichlet null
    VectorInterface h("h");
    ReInitVector(h);
    
    VectorInterface w("w");
    ReInitVector(w);

    GetSolver()->Equ(h,1.0,f);
    GetSolver()->AdjointForm(h,z,-1.0);
    GetSolver()->SetBoundaryVectorZero(h);
    
    Zero(w);
    info.reset();
    NewtonLinearSolve(w,h,info);
    GetSolver()->Add(z,1.0,w);
    
    DeleteVector(h);
    DeleteVector(w);
  }  

  //////////////////////////////////////////////////

  void ThetaMultiLevelSolver::SolveAdjointFirst(VectorInterface& z,  VectorInterface& f, CGInfo& info)
  {
    ThetaSolver* TS = dynamic_cast< ThetaSolver*> (GetSolver());
    assert(TS);
    
    Zero(z);
    VectorInterface h("h");
    ReInitVector(h);
    VectorInterface w("w");
    ReInitVector(w);
    
    GetSolver()->Equ(h,1.0,f);
    TS->GetMassMatrix()->vmult_time(TS->GetGV(h), TS->GetGV(z), TS->GetTimePattern(),-1.0);
    GetSolver()->SetBoundaryVectorZero(h);
    
    Zero(w);
    
    info.reset();
    NewtonLinearSolve(w,h,info);
    GetSolver()->Add(z,1.0,w);
    
    
    DeleteVector(h);
    DeleteVector(w);
  }  

  void ThetaSolver::RhsOld(VectorInterface& f, const VectorInterface& old,double time, double dt, double theta) 
  {
    SetTimeData(dt,theta,time,0,0);
  
    GetMassMatrix()->vmult_time(GetGV(f),GetGV(old), GetTimePattern(), 1.0);
    
    StdSolver::Rhs(f,(1.0 - theta) * dt);
    StdSolver::Form(f,old, - (1.0 - theta) * dt);
  }

  void ThetaSolver::RhsNew(VectorInterface& f, double time, double dt, double theta) 
  {
    SetTimeData(dt,theta,time,0,0);  
  
    StdSolver::Rhs(f, theta * dt);
  }


  /////////////

  void ThetaSolver::AdjointRhsOld(VectorInterface& f, const VectorInterface& zold, const VectorInterface& uold,double time, double dt, double theta)  
  {
    SetTimeData(dt,theta,time,0,0);
    
    GetMassMatrix()->vmult_time(GetGV(f),GetGV(zold), GetTimePattern(), 1.0);
    AddNodeVector("u",uold);
    StdSolver::AdjointForm(f,zold, - (1.0 - theta) * dt);
    DeleteNodeVector("u");
  }

  void ThetaSolver::AdjointRhsNew(VectorInterface& f, double time, double dt, double theta) 
  {
  }


  void ThetaSolver::AdjointRhsFunctional(VectorInterface& f, const GlobalVector& U, double s)
  {
    MyFuncRhs rhs_end_func;

    GetDiscretization()->AddNodeVector("u",&U);
    HNAverageData();
    
    GetDiscretization()->Rhs(GetGV(f),rhs_end_func,s);
    //GetDiscretization()->DiracRhs(GetGV(f),rhs_end_func,s);

    GetDiscretization()->DeleteNodeVector("u");
  }

  /*-------------------------------------------------------*/
  
  void ThetaSolver::RhsEndFunctional(VectorInterface& f, GlobalVector& U)
  {
    MyFuncRhs rhs_end_func;
    
    GetDiscretization()->AddNodeVector("u",&U);
    HNAverageData();
    
    
    GetDiscretization()->Rhs(GetGV(f),rhs_end_func,1.0);
    //GetDiscretization()->DiracRhs(GetGV(f),rhs_end_func,1.0);

    GetDiscretization()->DeleteNodeVector("u");
  }

  
  /*-------------------------------------------------------*/

  void ThetaSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
  {
    StdSolver::Form(gy,gx,d * _dt * _theta);

    const GlobalVector& x = GetGV(gx);
    GlobalVector& y = GetGV(gy);
    
    assert(y.n()==x.n());
    GetMassMatrix()->vmult_time(y,x,GetTimePattern(),d);
  }
  
  /*-------------------------------------------------------*/

  void ThetaSolver::FormOnly(VectorInterface& gy, const VectorInterface& gx, double d) const
  {
    StdSolver::Form(gy,gx,d);
  }

  /* ---------------------------------------- */

  void ThetaSolver::AdjointForm(VectorInterface& gy, const VectorInterface& gx, double d) const
  {
    StdSolver::AdjointForm(gy,gx,d * _dt * _theta);

    const GlobalVector& x = GetGV(gx);
    GlobalVector& y = GetGV(gy);
    
    assert(y.n()==x.n());
    GetMassMatrix()->vmult_time(y,x,GetTimePattern(),d);
  }

  /* -------------------------------------------------- */

  void ThetaSolver::AdjointFormOnly(VectorInterface& gy, const VectorInterface& gx, double d) const
  {
    StdSolver::AdjointForm(gy,gx,d);
  }

  /*-------------------------------------------------------*/
  
  void ThetaSolver::AssembleMatrix(const VectorInterface& gu, double d)
  {
    StdSolver::AssembleMatrix(gu,d * _dt * _theta);
    
    GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),GetTimePattern(),d);
    
    StdSolver::PeriodicMatrix();
    StdSolver::DirichletMatrix();
  }

  // --------------------------------------------------
  
  void ThetaSolver::AssembleAdjointMatrix(const VectorInterface& u, double dt, double theta) 
  {
    StdSolver::MatrixZero();
    
    StdSolver::AssembleMatrix(u,dt * theta);
    GetMatrix()->transpose();
        
    GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),GetTimePattern(),1.0);
    
    StdSolver::PeriodicMatrix();
    StdSolver::DirichletMatrix();

    ComputeIlu(u);
  }

  // --------------------------------------------------

  void ThetaSolver::AssembleAdjointMatrixFirst(const VectorInterface& u)
  {
    StdSolver::MatrixZero();
    
    GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),GetTimePattern(),1.0);
    
    StdSolver::PeriodicMatrix();
    StdSolver::DirichletMatrix();

    ComputeIlu(u);
  }
  
  
}
