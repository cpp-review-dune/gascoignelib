#ifndef __SolverInfos_h
#define __SolverInfos_h

#include "nlinfo.h"
#include "paramfile.h"
#include <map>

/*---------------------------------------------------------------*/

namespace Gascoigne
{
class SolverInfos
{
  protected:
    
    std::map<std::string,CGInfo*> __L;
    std::map<std::string,NLInfo*> _NL;

    std::string _linearsolve;

  public:
    SolverInfos() { }
    virtual ~SolverInfos();

    virtual void BasicInit(const ParamFile *param);
    CGInfo& GetLInfo(std::string s = "State") const;
    NLInfo& GetNLInfo(std::string s = "State") const;
};

/*---------------------------------------------------------------*/
}

#endif
