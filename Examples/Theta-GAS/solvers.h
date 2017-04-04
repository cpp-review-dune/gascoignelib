/*----------------------------   solvers.h     ---------------------------*/
/*      $Id: solvers.h,v 1.1 2012/09/02 16:43:20 richter Exp $                 */
#ifndef __solvers_H
#define __solvers_H
/*----------------------------   solvers.h     ---------------------------*/


#include "stdmultilevelsolver.h"
#include "stdtimesolver.h"


namespace Gascoigne
{

  class ThetaSolver : public StdTimeSolver
  {

  public:
    void TimeRhsOperator(VectorInterface& f, const VectorInterface& u) const {  abort();}
    void TimeRhs(int k, VectorInterface& f) const { abort(); }

    void RhsOld(VectorInterface& f, const VectorInterface& old,double time, double dt, double theta);
    void RhsNew(VectorInterface& f, double time, double dt, double theta);

    void AdjointRhsOld(VectorInterface& f, const VectorInterface& zold, const VectorInterface& uold,double time, double dt, double theta);
    void AdjointRhsNew(VectorInterface& f, double time, double dt, double theta);
	
    void RhsEndFunctional(VectorInterface& f, GlobalVector& U);
    void AdjointRhsFunctional(VectorInterface& f, const GlobalVector& U, double s);
    
    void AssembleMatrix(const VectorInterface& gu, double d);

    void AssembleAdjointMatrix(const VectorInterface& u, double dt, double theta) ;
    void AssembleAdjointMatrixFirst(const VectorInterface& u);
    void FormOnly(VectorInterface& gy, const VectorInterface& gx, double d) const    ;
    void AdjointFormOnly(VectorInterface& gy, const VectorInterface& gx, double d) const;
    void Form(VectorInterface& gy, const VectorInterface& gx, double d) const;
    void AdjointForm(VectorInterface& gy, const VectorInterface& gx, double d) const;

    
  };
  
  class ThetaMultiLevelSolver : public StdMultiLevelSolver
  {
  public:
    SolverInterface* NewSolver(int solverlevel)
    { return new ThetaSolver; }

    void SolveAdjoint(VectorInterface& z, VectorInterface& f, CGInfo& info);
    void SolveAdjointFirst(VectorInterface& z, VectorInterface& f, CGInfo& info);
  };
  
}


/*----------------------------   solvers.h     ---------------------------*/
/* end of #ifndef __solvers_H */
#endif
/*----------------------------   solvers.h     ---------------------------*/
