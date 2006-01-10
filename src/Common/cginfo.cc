#include  "cginfo.h"
#include  "math.h"
#include  "gascoignemath.h"
#include  <iostream>

using namespace std;

/*******************************************************************/

namespace Gascoigne
{
ostream& operator<<(ostream &s, const StatisticData& A)
{
  s << "StatisticData\n";
  s << "rate    " << "\t"<< A.rate() << endl;
  s << "lastrate" <<"\t"<< A.lastrate() << endl;
  s << "totaliter" <<"\t"<< A.totaliter() << endl;
  return s;
}

/*******************************************************************/

ostream& operator<<(ostream &s, const ControlData& A)
{
  s << "ControlData\n";
  s << "status" <<"\t"<< A.status() << endl;
  s << "firstresidual" <<"\t"<< A.firstresidual() << endl;
  s << "residual" <<"\t"<< A.residual() << endl;
  s << "aimedresidual" <<"\t"<< A.aimedresidual() << endl;
  s << "previousresidual" << A.previousresidual() << endl;
  s << "iteration" <<"\t"<< A.iteration() << endl;
  s << "correction" <<"\t"<< A.correction() << endl;
  s << "previouscorrection" <<"\t"<< A.previouscorrection() << endl;
  return s;
}

/*******************************************************************/

ostream& operator<<(ostream &s, const UserData& A)
{
  s << "UserData\n";
  s << "tol" <<"\t"<< A.tol() << endl;
  s << "globaltol" <<"\t"<< A.globaltol() << endl;
  s << "breaktol" <<"\t"<< A.breaktol() << endl;
  s << "miniter" <<"\t"<< A.miniter() << endl;
  s << "maxiter" <<"\t"<< A.maxiter() << endl;
  s << "printstep" <<"\t"<< A.printstep() << endl;
  return s;
}

/*******************************************************************/

ostream& operator<<(ostream &s, const CGInfo& A)
{
  s << "CGInfo\n";
  s << A.statistics()<<endl;
  s << A.control()<<endl;
  s << A.user()<<endl;
  return s;
}

/*******************************************************************/

void StatisticData::reset()
{
  _totaliter = 0;
  _rate      = 0.;
  _lastrate  = 1.;
}

/*******************************************************************/

ControlData::ControlData()
{
  reset();
}

/*******************************************************************/

void ControlData::reset()
{
  _residual  = 0.;
  _correction  = 0.;
  _iteration = 0;
  _status    = "waiting";
  _previousresidual = 1.;
  _previouscorrection = 1.;
}

/*******************************************************************/

void CGInfo::compute_reduction_rate()
{
  double b = CD.residual()/CD.firstresidual();
  double p = 1./max_int(1,CD.iteration());
  SD.rate() = pow(b,p);  
}

/*******************************************************************/

CGInfo::CGInfo(double f, double t, int p, int m, const string& txt)
{
  UD.text() = txt;
  
  UD.miniter  () = 0;
  UD.maxiter  () = m;
  UD.globaltol() = t;
  UD.printstep() = p;
  UD.tol()       = f;
  UD.breaktol()  = 1.e15;
 
  CD.reset();
  //SD.reset();
}

/*******************************************************************/

CGInfo::CGInfo(const string& txt)
{
  UD.text() = txt;
  double f = 1.e-6;
  double t=1.e-14;
  int p = 10; 
  int min =   0;
  int max = 100;
  
  UD.miniter  () = min;
  UD.maxiter  () = max;
  UD.globaltol() = t;
  UD.printstep() = p;
  UD.tol()       = f;
  UD.breaktol()  = 1.e15;
 
  CD.reset();
  SD.reset();
}

/*------------------------------------------------------------------*/

void CGInfo::reset()
{
  CD.reset();
  SD.rate() = 0.;
  SD.lastrate() = 0.;
}

/*------------------------------------------------------------------*/

bool CGInfo::check(double resi, double cori)
{
  double res=fabs(resi);
  double cor=fabs(cori);

  CD.residual  () = res;
  CD.correction() = cor;

  if (CD.status()!="running")
    {
      CD.firstresidual()    = res;
      CD.previousresidual() = res;

      CD.previouscorrection() = cor;
      CD.aimedresidual() = res*UD.tol();
      CD.status()        = "running";
    }
  else
    {
      SD.totaliter()++;
      SD.lastrate() = res/CD.previousresidual();
      compute_reduction_rate();
      CD.previousresidual() = res;
    }
  if ( CD.residual()< max(CD.aimedresidual(),UD.globaltol()) )
    {
      CD.status() = "converged";
    }
  else if ( CD.iteration()>=UD.maxiter() )
    {
      CD.status() = "too many iterations";
    }
  else if ( CD.residual()>UD.breaktol() )
    {
      CD.status() = "exploded";
    }
  else if ( CD.residual()>CD.previousresidual() )
    {
      CD.status() = "diverging";
    }
  if (UD.printstep() && !(CD.iteration()%UD.printstep()) )
    {
      int prec = cout.precision();
      cout.precision(5);
      cout << UD.text() << " " << CD.iteration() << " Res: " << CD.residual();
      cout << " Cor: " << CD.correction()<< endl;
      cout.precision(prec);
    }
  CD.iteration()++;
  if ( CD.status()    == "running"    )  return 0;
  if ( CD.iteration() <  UD.miniter() )  return 0;

  return 1;
}
}
