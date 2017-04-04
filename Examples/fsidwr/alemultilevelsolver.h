/*----------------------------   alemultilevelsolver.h     ---------------------------*/
/*      $Id: alemultilevelsolver.h,v 1.2 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __alemultilevelsolver_H
#define __alemultilevelsolver_H
/*----------------------------   alemultilevelsolver.h     ---------------------------*/

#include "stdmultilevelsolver.h"
#include "alesolver.h"

namespace Gascoigne
{
  
  class AleMultiLevelSolver : public StdMultiLevelSolver
  {
  public:
    std::string GetName() const {return "ALE MultiLevelSolver";}

    SolverInterface* NewSolver(int solverlevel)
    { return new AleSolver; }

    const AleSolver* GetAleSolver(int l) const
    {
      assert(dynamic_cast<const AleSolver*> (GetSolver(l)));
      return dynamic_cast<const AleSolver*> (GetSolver(l));
    }
    
    void NewMgInterpolator();
    
    void NewtonLinearSolve(VectorInterface& x, const VectorInterface& b, CGInfo& info);
    
    void smooth(int l, int iter, VectorInterface& x, VectorInterface& b, VectorInterface& h);
    void mgstep(vector<double>& res, vector<double>& rw, 
		int l, int finelevel, int coarselevel, string& p0, string p,
		VectorInterface& u, VectorInterface& b, VectorInterface& v);
    

/*     virtual void LinearMg_f(int minlevel, int maxlevel, VectorInterface& u, const VectorInterface& f, CGInfo&); */
/*     virtual void LinearMg_s(int minlevel, int maxlevel, VectorInterface& u, const VectorInterface& f, CGInfo&); */

    
/*     virtual void mgstep_f(std::vector<double>& res, std::vector<double>& rw, int l, int maxl, int minl, std::string& p0, std::string p, VectorInterface& u, VectorInterface& b, VectorInterface& v); */
/*     virtual void mgstep_s(std::vector<double>& res, std::vector<double>& rw, int l, int maxl, int minl, std::string& p0, std::string p, VectorInterface& u, VectorInterface& b, VectorInterface& v); */

  };

}



/*----------------------------   alemultilevelsolver.h     ---------------------------*/
/* end of #ifndef __alemultilevelsolver_H */
#endif
/*----------------------------   alemultilevelsolver.h     ---------------------------*/
