#include "multilevelsolverdata.h"
#include "filescanner.h"
#include "stringutil.h"

using namespace std;

/**********************************************************/

namespace Gascoigne
{
MultiLevelSolverData::~MultiLevelSolverData()
{
}

/**********************************************************/

void MultiLevelSolverData::BasicInit(const ParamFile *param)
{
  _countresidual = 0; 
  
  DataFormatHandler DFH;
  DFH.insert("linearsolve",         &_linearsolve,        "mg");
  DFH.insert("nonlinearsolve",      &_nonlinearsolve,     "newton");
  DFH.insert("solver",              &_solver,             "stat");
  DFH.insert("mgomega",             &_mgomega,            1.);
  DFH.insert("coarselevel",         &_coarselevel,        0);
  DFH.insert("mgtype",              &_mgtype,             "V");
  DFH.insert("projection",          &_projectionflag,     0);

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(param,"Multilevelsolver");

  if ((_mgtype!="V") && (_mgtype!="W") && (_mgtype!="F"))
  {
    _mgtype = "V";
  }

  vector<string> v = StringSplit(_linearsolve.c_str(),'_');
  if( (v.size()==2) && (v[0]=="gmres") )
  {
    _gmresmemsize = atoi(v[1].c_str());
    _linearsolve = "gmres";
  }
}
}
