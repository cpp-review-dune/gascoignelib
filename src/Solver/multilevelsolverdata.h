#ifndef  __multilevelsolverdata_h
#define  __multilevelsolverdata_h

#include  "nlinfo.h"

/*-------------------------------------------------------------*/

class MultiLevelSolverData
{
public:

  std::string    linearsolve, nonlinearsolve;
  std::string    mgtype, solver;
  int        coarselevel;
  double     mgomega;
  CGInfo     info, dualinfo, precinfo;
  NLInfo     nlinfo, dualnlinfo;

  int        countresidual;
  int        gmresmemsize, projectionflag;

  // constructors

  MultiLevelSolverData(const std::string& filename);
  ~MultiLevelSolverData();
};

#endif
