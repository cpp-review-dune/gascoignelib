#ifndef __MultiLevelSolverData_h
#define __MultiLevelSolverData_h

#include "nlinfo.h"
#include "cginfo.h"
#include "paramfile.h"
#include <map>
#include <iostream>

/**********************************************************/

class MultiLevelSolverData
{
  private:

  protected:
    std::string _solver, _mgtype, _linearsolve, _nonlinearsolve;
    int _countresidual, _coarselevel, _gmresmemsize, _projectionflag;
    double _mgomega;
    
    std::map<std::string,CGInfo *> _L;
    std::map<std::string,NLInfo *> _NL;

  public:
    MultiLevelSolverData() { }
    virtual ~MultiLevelSolverData();

    virtual void BasicInit(const Gascoigne::ParamFile *param);
    CGInfo &GetLInfo(std::string s = "State") const;
    NLInfo &GetNLInfo(std::string s = "State") const;

    std::string &Solver() { return _solver; }
    int &CountResidual() { return _countresidual; }
    int &CoarseLevel() { return _coarselevel; }
    std::string &MgType() { return _mgtype; }
    double &MgOmega() { return _mgomega; }
    std::string &LinearSolve() { return _linearsolve; }
    std::string &NonLinearSolve() { return _nonlinearsolve; }
    int &GmresMemSize() { return _gmresmemsize; }
    int &ProjectionFlag() { return _projectionflag; }

};

/**********************************************************/

inline CGInfo &MultiLevelSolverData::GetLInfo(std::string s) const
{
  std::map<std::string,CGInfo *>::const_iterator iter = _L.find(s);
  if(iter != _L.end())
  {
    return *iter->second; 
  }
  else
  {
    std::cerr << "No such LinearInfo: \"" << s << "\"!" << std::endl;
    abort();
  }
}

inline NLInfo &MultiLevelSolverData::GetNLInfo(std::string s) const
{
  std::map<std::string,NLInfo *>::const_iterator iter = _NL.find(s);
  if(iter != _NL.end())
  {
    return *iter->second; 
  }
  else
  {
    std::cerr << "No such NonLinearInfo: \"" << s << "\"!" << std::endl;
    abort();
  }
}

/**********************************************************/

#endif
