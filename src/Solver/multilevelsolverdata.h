#ifndef __MultiLevelSolverData_h
#define __MultiLevelSolverData_h

#include "paramfile.h"
#include <map>
#include <iostream>

/**********************************************************/

namespace Gascoigne
{
class MultiLevelSolverData
{
  protected:

    std::string  _solver, _mgtype, _linearsolve, _nonlinearsolve;
    int          _countresidual, _coarselevel, _gmresmemsize, _projectionflag;
    double       _mgomega;
    
  public:

    MultiLevelSolverData() { }
    virtual ~MultiLevelSolverData();

    virtual void BasicInit(const ParamFile *param);

    std::string& Solver()        { return _solver; }
    int&         CountResidual() { return _countresidual; }
    int&         CoarseLevel()   { return _coarselevel; }
    std::string& MgType()        { return _mgtype; }
    double&      MgOmega()       { return _mgomega; }
    std::string& LinearSolve()   { return _linearsolve; }
    std::string& NonLinearSolve(){ return _nonlinearsolve; }
    int&         GmresMemSize()  { return _gmresmemsize; }
    int&         ProjectionFlag(){ return _projectionflag; }
};

/**********************************************************/

}

/**********************************************************/

#endif
