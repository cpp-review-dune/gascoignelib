#include "solverinfos.h"
#include <cassert>
#include "filescanner.h"

using namespace std;

namespace Gascoigne
{

SolverInfos::~SolverInfos()
{
  map<string,CGInfo *>::iterator p = __L.begin();
  while(p != __L.end())
  { 
    if(p->second) 
    {
      delete p->second; 
      p->second = NULL;
    } 
    p++;
  }
  
  map<string,NLInfo *>::iterator q = _NL.begin();
  while(q != _NL.end())
  { 
    if(q->second) 
    {
      delete q->second; 
      q->second = NULL;
    } 
    q++;
  }
}

/*---------------------------------------------------------------*/

CGInfo& SolverInfos::GetLInfo(std::string s) const
{
  map<string,CGInfo *>::const_iterator iter = __L.find(s);
  assert(iter!=__L.end());
  return *iter->second; 
}

/*---------------------------------------------------------------*/

NLInfo& SolverInfos::GetNLInfo(std::string s) const
{
  map<string,NLInfo *>::const_iterator iter = _NL.find(s);
  assert(iter!=_NL.end());
  return *iter->second; 
}

/*---------------------------------------------------------------*/

void SolverInfos::BasicInit(const ParamFile *param)
{
  __L["State"]   = new CGInfo();
  _NL["State"]   = new NLInfo(GetLInfo());
  __L["Precond"] = new CGInfo();
  
  DataFormatHandler DFH;

  double linear_tol, linear_globaltol;
  int    linear_maxiter, linear_pstep;

  double nonlinear_tol, nonlinear_globaltol, nonlinear_rho, nonlinear_increase;
  int    nonlinear_maxiter, nonlinear_pstep, nonlinear_damp;

  double prec_tol, prec_globaltol;
  int    prec_maxiter, prec_pstep;

  string linearsolve;

  DFH.insert("linearsolve",         &linearsolve,        "mg");
  DFH.insert("linear_tol",          &linear_tol,          1.e-2);
  DFH.insert("linear_globaltol",    &linear_globaltol,    1.e-12);
  DFH.insert("linear_maxiter",      &linear_maxiter,      10);
  DFH.insert("linear_pstep",        &linear_pstep,        0);

  DFH.insert("nonlinear_tol",       &nonlinear_tol,       1.e-4);
  DFH.insert("nonlinear_globaltol", &nonlinear_globaltol, 1.e-12);
  DFH.insert("nonlinear_rho",       &nonlinear_rho,       0.001);
  DFH.insert("nonlinear_maxiter",   &nonlinear_maxiter,   10);
  DFH.insert("nonlinear_pstep",     &nonlinear_pstep,     0);
  DFH.insert("nonlinear_damp",      &nonlinear_damp,      4);
  DFH.insert("nonlinear_increase",  &nonlinear_increase,  1.e3);

  DFH.insert("prec_tol",            &prec_tol,            1.e-12);
  DFH.insert("prec_globaltol",      &prec_globaltol,      1.e-12);
  DFH.insert("prec_maxiter",        &prec_maxiter,        1);
  DFH.insert("prec_pstep",          &prec_pstep,          0);

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(param,"Multilevelsolver");

  GetLInfo().user().tol()                = linear_tol;
  GetLInfo().user().globaltol()          = linear_globaltol;
  GetLInfo().user().maxiter  ()          = linear_maxiter;
  GetLInfo().user().printstep()          = linear_pstep;
  if (linearsolve=="mg")
    {
      GetLInfo().user().text() = "MGInfo";
    }
  else
    {
      GetLInfo().user().text() = "CGInfo";
    }
                                         
  GetNLInfo().user().tol()               = nonlinear_tol;
  GetNLInfo().user().globaltol()         = nonlinear_globaltol;
  GetNLInfo().user().rho()               = nonlinear_rho;
  GetNLInfo().user().maxiter()           = nonlinear_maxiter;
  GetNLInfo().user().printstep()         = nonlinear_pstep;
  GetNLInfo().user().maxrelax()          = nonlinear_damp;
  GetNLInfo().user().maxresincrease()    = nonlinear_increase;

  GetLInfo("Precond").user().tol()       = prec_tol;
  GetLInfo("Precond").user().globaltol() = prec_globaltol;
  GetLInfo("Precond").user().maxiter()   = prec_maxiter;
  GetLInfo("Precond").user().printstep() = prec_pstep;
  GetLInfo("Precond").user().text()      = "PrecInfo";
}

/*---------------------------------------------------------------*/
}

