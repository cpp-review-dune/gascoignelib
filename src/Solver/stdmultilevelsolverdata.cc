#include "stdmultilevelsolverdata.h"
#include "filescanner.h"

using namespace std;

/**********************************************************/

namespace Gascoigne
{
StdMultiLevelSolverData::~StdMultiLevelSolverData()
{
}

/**********************************************************/

void StdMultiLevelSolverData::BasicInit(const ParamFile *param)
{
  _countresidual = 0; 
  
  DataFormatHandler DFH;
  DFH.insert("linearsolve",         &_linearsolve,        "mg");
  DFH.insert("nonlinearsolve",      &_nonlinearsolve,     "newton");
  DFH.insert("solver",              &_solver,             "stat");
  DFH.insert("mgomega",             &_mgomega,            1.);
  DFH.insert("coarselevel",         &_coarselevel,        0);
  DFH.insert("mgtype",              &_mgtype,             "V");

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(param,"Multilevelsolver");

  if ((_mgtype!="V") && (_mgtype!="W") && (_mgtype!="F"))
  {
    _mgtype = "V";
  }
}
}
