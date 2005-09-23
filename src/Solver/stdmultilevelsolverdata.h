#ifndef __StdMultiLevelSolverData_h
#define __StdMultiLevelSolverData_h

#include "paramfile.h"
#include "cginfo.h"
#include <map>
#include <iostream>

/**********************************************************/

namespace Gascoigne
{
class StdMultiLevelSolverData
{
  protected:

    std::string  _solver, _mgtype, _linearsolve, _nonlinearsolve;
    int          _countresidual, _coarselevel;
    int          _i_show_nonlinear_comp_residuals,_i_show_linear_comp_residuals,_i_show_comp_residual_names;
    int          _i_save_nonlinear_comp_residuals,_i_save_linear_comp_residuals;
    double       _mgomega;
    
    int        _gmresmemsize;
    CGInfo     precinfo;
  
  public:

    StdMultiLevelSolverData() { }
    virtual ~StdMultiLevelSolverData();

    virtual void BasicInit(const ParamFile *param);

    std::string& Solver()        { return _solver; }
    int&         CountResidual() { return _countresidual; }
    int&         CoarseLevel()   { return _coarselevel; }
    std::string& MgType()        { return _mgtype; }
    double&      MgOmega()       { return _mgomega; }
    std::string& LinearSolve()   { return _linearsolve; }
    std::string& NonLinearSolve(){ return _nonlinearsolve; }

    int  GmresMemSize()    const { return _gmresmemsize; }
          CGInfo& GetPrecInfo()       { return precinfo;}
    const CGInfo& GetPrecInfo() const { return precinfo;}

    int  SaveNonLinearCompResiduals() const { return _i_save_nonlinear_comp_residuals; }
    int  SaveLinearCompResiduals()    const { return _i_save_linear_comp_residuals;    }
    int  ShowNonLinearCompResiduals() const { return _i_show_nonlinear_comp_residuals; }
    int  ShowLinearCompResiduals()    const { return _i_show_linear_comp_residuals;    }
    int  ShowCompResidualNames()      const { 
      return ( _i_show_comp_residual_names && (_i_show_nonlinear_comp_residuals || _i_show_linear_comp_residuals ));
    }
};

/**********************************************************/

}

/**********************************************************/

#endif
