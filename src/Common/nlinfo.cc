#include "nlinfo.h"
#include "math.h"
#include "gascoignemath.h"
#include  <iostream>

using namespace std;

namespace Gascoigne
{

/*******************************************************************/

ostream& operator<<(ostream &s, const NLStatisticData& A)
{
  s << "NLStatisticData\n";
  s << "newmatrix" <<"\t"<< A.newmatrix()<< endl;
  s << "totalmatrix" <<"\t"<< A.totalmatrix()<< endl;
  s << static_cast<StatisticData>(A);
  return s;
}

/*******************************************************************/

ostream& operator<<(ostream &s, const NLControlData& A)
{
  s << "NLControlData\n";
  s << "relax" <<"\t"<< A.relax()<< endl;
  s << "newmatrix" <<"\t"<< A.newmatrix()<< endl;
  s << "laststepbad" <<"\t"<< A.laststepbad()<< endl;
  s << static_cast<ControlData>(A);
  return s;
}

/*******************************************************************/

ostream& operator<<(ostream &s, const NLUserData& A)
{
  s << "NLUserData\n";
  s << "rho" <<"\t"<< A.rho()<< endl;
  s << "linrho" <<"\t"<< A.linrho()<< endl;
  s << "maxrelax" <<"\t"<< A.maxrelax()<< endl;
  s << static_cast<UserData>(A);
  return s;
}

/*******************************************************************/

ostream& operator<<(ostream &s, const NLInfo& A)
{
  s << "NLInfo\n";
  s << A.statistics()<<endl;
  s << A.control()<<endl;
  s << A.user()<<endl;
  return s;
}

/*******************************************************************/

void NLStatisticData::reset()
{
  StatisticData::reset();
  _totalmatrix = 0;
  _newmatrix = 0;
}

/*******************************************************************/

NLControlData::NLControlData()
{
  _matrixmustbebuild = 1;
}

/*******************************************************************/

void NLControlData::reset()
{
  ControlData::reset();
  _relax     = 0;
  _laststepbad = 0;
//   _matrixmustbebuild=1;
  //_newmatrix = 1;
}

/*******************************************************************/

NLUserData::NLUserData() : UserData()
{
  breaktol()  = 1.e15;
  maxrelax()  = 11;
  rho()       = 0.3;
  linrho()    = 1.;
  maxresincrease() = 1.e3;
}

/*******************************************************************/

void NLInfo::new_matrix()
{
  CD.newmatrix() = 1;
}

/*******************************************************************/

NLInfo::NLInfo(CGInfo& info, double f, double t, int p, int m, const string& txt) :
Linfo(info)
{
  UD.text() = txt;
  
  UD.maxiter  () = m;
  UD.globaltol() = t;
  UD.printstep() = p;
  UD.tol()       = f;
 
  CD.reset();
  SD.reset();
}

/*******************************************************************/

void NLInfo::reset()
{
  SD.newmatrix() = 0;
  CD.reset();
  //SD.reset();
  Linfo.reset();
  SD.rate() = 0.;
}

/*******************************************************************/

void NLInfo::compute_reduction_rate()
{
  double b = CD.residual()/CD.firstresidual();
  double p = 1./max_int(1,CD.iteration());
  SD.rate() = pow(b,p);  

  if (CD.iteration()>1)
    {
      SD.lastrate() = CD.residual()/CD.previousresidual();
    }
  else
    SD.lastrate() = SD.rate();
}

/*******************************************************************/
 
void NLInfo::matrix_control()
{
//   cerr << "NLInfo::matrix_control()\t";
//   cerr << GetLinearInfo().statistics().lastrate()<<"\t"<<UD.linrho()<<endl;

  bool nonlinbad = (SD.lastrate()>UD.rho());
  bool linbad    = (GetLinearInfo().statistics().lastrate()>UD.linrho());

  CD.newmatrix() = 0;
  if( nonlinbad || linbad )
    {
      CD.newmatrix() = 1;
      SD.totalmatrix()++;
      SD.newmatrix()++;
    }
}

/*******************************************************************/
 
string NLInfo::check_damping(int dampit, double res)
{
  CD.residual() = fabs(res);
  CD.relax()    = dampit;

  if (CD.residual()>1.e20)     return "exploded"; 
  if (!(CD.residual()<=1.e20)) return "exploded"; 

  if (CD.residual()<=CD.previousresidual()) return "ok";
  return "continue";
}

/*******************************************************************/
 
bool NLInfo::check(double resi, double cori)
{
  return check(CD.iteration(),resi,cori);
}

/*******************************************************************/
 
bool NLInfo::check(int iter, double resi, double cori)
{
  double res=fabs(resi);

  double cor=fabs(cori);
  bool newiteration = 0;

  if (CD.status()=="diverged")   return 1;
  if (CD.status()=="stagnation") return 1;

  if (CD.status()=="waiting")
    {
      newiteration = 1;
    }
  else if (iter>CD.iteration())
    {
      newiteration = 1;
      CD.iteration()++;
      SD.totaliter()++;
    }

  CD.status() = "running";

  int thisstepbad = 0;

  if (!CD.iteration())
    {
      CD.residual () = res;
      CD.previousresidual() = res;
      CD.correction() = cor;
      CD.firstresidual() = res;
      CD.aimedresidual() = res*UD.tol();
    }
  else
    {
      double r0 = CD.residual();
      double c0 = CD.correction();
      CD.residual () = res;
      CD.correction() = cor;
      compute_reduction_rate();
      if (CD.residual()>=2.*CD.previousresidual())
	{
	  thisstepbad = 1;
	}
      if (newiteration) 
	{
	  CD.previousresidual()   = r0;
	  CD.previouscorrection() = c0;
	}
    }

  if (thisstepbad && CD.laststepbad())
    {
      CD.status() = "diverged";
    }
  else if ( CD.residual()< max(CD.aimedresidual(),UD.globaltol()) )
    {
      CD.status() = "converged";
    }
  else if ( (CD.iteration()>=UD.maxiter()) && (SD.rate()<1.) )//<0.95) )
    {
      CD.status() = "slow_convergence";
    }
  else if ( CD.iteration()>=UD.maxiter() )
    {
      CD.status() = "stagnation";
    }
  else if ( CD.residual()>UD.breaktol() )
    {
      CD.status() = "exploded";
    }
  else if (CD.iteration() && SD.lastrate()> UD.maxresincrease())
    {
      CD.status() = "exploded";
    }
  if (UD.printstep() && !(CD.iteration()%UD.printstep()) )
    {
      cout.setf(ios::scientific,ios::floatfield);
      cout.precision(5);
      cout << UD.text() << " " << CD.iteration() << "\t" << CD.residual();
      cout << " [" << CD.correction()<< "] ";
      cout << Linfo.control().residual();
      cout << " [" << Linfo.control().correction() << "] ";
      cout << "\t" << Linfo.control().iteration() << "\t";
      cout << Linfo.statistics().rate() << endl;
    }

  CD.laststepbad() = thisstepbad;
  matrix_control();

  if (CD.status()=="running") return 0;

  SD.totaliter() += CD.iteration();

  return 1;
}
}
