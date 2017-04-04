/*----------------------------   loop.h     ---------------------------*/
/*      $Id: loop.h,v 1.3 2012/09/19 14:43:03 richter Exp $                 */
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/


#include "stdtimeloop.h"





namespace Gascoigne
{

  
  class ThetaLoop : public StdTimeLoop
  {
    double __Tstart,__Tstop, __time, __theta0;
    nvector<double> __DT, __THETA, __TIME;
    nvector<GlobalVector> __primal_solution;
    nvector<GlobalVector> __dual_solution;
    int reflevel;

  public:
    double EstimateInt   (DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f);
    double EstimatePrimal(DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f);
    double EstimateDual  (DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f);
    double EstimatePrimalTheta(DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f);
    double EstimateDualTheta  (DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f);
    double Estimator(double& ei, double& ep, double& ed, DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f);
    double EstimatorTheta(double& ei, double& ep, double& ed, DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f);
    
    double GetEndFunctional(VectorInterface& f,VectorInterface& u);
    double GetTimeFunctional(VectorInterface& f,VectorInterface& u); 
    
    void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC=NULL);
    void PrimalLoop(VectorInterface& u, VectorInterface& f);
    void AdjointLoop(VectorInterface& z, VectorInterface& f, VectorInterface& u);

   void SetTimeData(double time, double dt, double theta);
   void run(const std::string& problemlabel, int NITER, double THETA);
   void run_theta(const std::string& problemlabel, int NITER, double THETA);
   void run_adaptive_theta(const std::string& problemlabel, int NITER, double THETA);
   void run_adaptive(const std::string& problemlabel, int NITER, double THETA);
   void run_adaptive_theta2(const std::string& problemlabel, int NITER, double THETA);
   void run_exact(const std::string& problemlabel, int NITER, double THETA);
   void Extra(std::vector<double>& f);
   void SetReflevel(int r)
   {
     reflevel=r;
    }

   
  };

}


/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
