#include  "multilevelsolverdata.h"
#include  "filescanner.h"
#include  "stringutil.h"


using namespace std;
using namespace Gascoigne;

/*-------------------------------------------------------------*/

MultiLevelSolverData::MultiLevelSolverData(const ParamFile* pf) 
  : info("Lin"), dualinfo("Dual"), nlinfo(info), dualnlinfo(dualinfo), 
    countresidual(0)
{
  DataFormatHandler DFH;
  DFH.insert("linearsolve"   ,&linearsolve   ,"mg");
  DFH.insert("nonlinearsolve",&nonlinearsolve,"newton");
  DFH.insert("solver"    ,&solver    ,"stat");
  DFH.insert("mgomega"       ,&mgomega    ,1.);

  double linear_tol, linear_globaltol;
  int    linear_maxiter, linear_pstep;

  double nonlinear_tol, nonlinear_globaltol,nonlinear_rho;
  int    nonlinear_maxiter, nonlinear_pstep;

  double dual_tol, dual_globaltol;
  int    dual_maxiter, dual_pstep;

  double prec_tol, prec_globaltol;
  int    prec_maxiter, prec_pstep;

  DFH.insert("coarselevel"      ,&coarselevel,0);
  
  DFH.insert("linear_tol"       ,&linear_tol,1.e-2);
  DFH.insert("linear_globaltol" ,&linear_globaltol,1.e-12);
  DFH.insert("linear_maxiter"   ,&linear_maxiter,10);
  DFH.insert("linear_pstep"     ,&linear_pstep,0);

  DFH.insert("nonlinear_tol"      ,&nonlinear_tol,1.e-4);
  DFH.insert("nonlinear_globaltol",&nonlinear_globaltol,1.e-12);
  DFH.insert("nonlinear_rho"      ,&nonlinear_rho,0.001);
  DFH.insert("nonlinear_maxiter"  ,&nonlinear_maxiter,10);
  DFH.insert("nonlinear_pstep"    ,&nonlinear_pstep,0);
  DFH.insert("nonlinear_damp"     ,&nlinfo.user().maxrelax(),4);
  DFH.insert("nonlinear_increase" ,&nlinfo.user().maxresincrease());

  DFH.insert("nonlinear_dualtol",&dualnlinfo.user().tol(),1.e-2);
  DFH.insert("dual_tol"       ,&dual_tol,1.e-2);
  DFH.insert("dual_globaltol" ,&dual_globaltol,1.e-12);
  DFH.insert("dual_maxiter"   ,&dual_maxiter,5);
  DFH.insert("dual_pstep"     ,&dual_pstep,0);

  DFH.insert("prec_tol"       ,&prec_tol,1.e-12);
  DFH.insert("prec_globaltol" ,&prec_globaltol,1.e-12);
  DFH.insert("prec_maxiter"   ,&prec_maxiter,1);
  DFH.insert("prec_pstep"     ,&prec_pstep,0);

  DFH.insert("mgtype"       ,&mgtype  ,"V");
  DFH.insert("projection"   ,&projectionflag,0);

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(pf,"Multilevelsolver");

  if ((mgtype!="V") && (mgtype!="W") && (mgtype!="F"))
    {
      mgtype = "V";
    }

  {
    vector<string> v = StringSplit(linearsolve.c_str(),'_');
    if( (v.size()==2) && (v[0]=="gmres") )
      {
	gmresmemsize = atoi(v[1].c_str());
	linearsolve = "gmres";
      }
  }

  info  .user().tol()       = linear_tol;
  info  .user().globaltol() = linear_globaltol;
  info  .user().maxiter  () = linear_maxiter;
  info  .user().printstep() = linear_pstep;

  nlinfo.user().tol()       = nonlinear_tol;
  nlinfo.user().globaltol() = nonlinear_globaltol;
  nlinfo.user().rho()       = nonlinear_rho;
  nlinfo.user().maxiter  () = nonlinear_maxiter;
  nlinfo.user().printstep() = nonlinear_pstep;

  dualinfo  .user().tol()       = dual_tol;
  dualinfo  .user().globaltol() = dual_globaltol;
  dualinfo  .user().maxiter  () = dual_maxiter;
  dualinfo  .user().printstep() = dual_pstep;

  dualnlinfo.user().tol()       = nonlinear_tol;
  //  dualnlinfo.user().tol()       = nonlinear_tol;
  dualnlinfo.user().globaltol() = nonlinear_globaltol;
  dualnlinfo.user().rho()       = 1.;
  dualnlinfo.user().maxiter  () = nonlinear_maxiter;
  dualnlinfo.user().printstep() = nonlinear_pstep;

  precinfo  .user().tol()       = prec_tol;
  precinfo  .user().globaltol() = prec_globaltol;
  precinfo  .user().maxiter  () = prec_maxiter;
  precinfo  .user().printstep() = prec_pstep;
  precinfo  .user().text()      = "PrecInfo";
  if (linearsolve=="mg")
    {
      info.user().text() = "MGInfo";
    }
}

MultiLevelSolverData::~MultiLevelSolverData()
{
}
