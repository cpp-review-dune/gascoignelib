#include "stdmultilevelsolverdata.h"
#include "filescanner.h"
#include "stringutil.h"

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
  
  double prec_tol, prec_globaltol;
  int    prec_maxiter, prec_pstep;

  DataFormatHandler DFH;
  DFH.insert("linearsolve",         &_linearsolve,        "mg");
  DFH.insert("nonlinearsolve",      &_nonlinearsolve,     "newton");
  DFH.insert("solver",              &_solver,             "stat");
  DFH.insert("mgomega",             &_mgomega,            1.);
  DFH.insert("coarselevel",         &_coarselevel,        0);
  DFH.insert("mgtype",              &_mgtype,             "V");

  DFH.insert("prec_tol",      &prec_tol,      1.e-12);
  DFH.insert("prec_globaltol",&prec_globaltol,1.e-12);
  DFH.insert("prec_maxiter",  &prec_maxiter,  1);
  DFH.insert("prec_pstep",    &prec_pstep,    0);
  DFH.insert("gmressize",     &_gmresmemsize, 10);

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(param,"Multilevelsolver");

  if ((_mgtype!="V") && (_mgtype!="W") && (_mgtype!="F"))
  {
    _mgtype = "V";
  }
  precinfo.user().tol()       = prec_tol;
  precinfo.user().globaltol() = prec_globaltol;
  precinfo.user().maxiter()   = prec_maxiter;
  precinfo.user().printstep() = prec_pstep;
  precinfo.user().text()      = "PrecInfo";

  std::vector<std::string> v = StringSplit(_linearsolve.c_str(),'_');
  if( (v.size()==2) && (v[0]=="gmres") )
  {
    _gmresmemsize = atoi(v[1].c_str());
    _linearsolve = "gmres";
  }
}
}
