#ifndef __StdMultiLevelSolverData_h
#define __StdMultiLevelSolverData_h

#include "paramfile.h"
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
    double       _mgomega;
    
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
};

/**********************************************************/

}

/**********************************************************/

#endif