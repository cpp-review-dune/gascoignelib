#include  "cginfo.h"
#include  "math.h"
#include  "fadamath.h"
#include  <iostream>

/*******************************************************************/

std::ostream& operator<<(std::ostream &s, const StatisticData& A)
{
  s << "StatisticData\n";
  s << "rate" << "\t"<< A.rate() << std::endl;
  s << "lastrate" <<"\t"<< A.lastrate() << std::endl;
  s << "totaliter" <<"\t"<< A.totaliter() << std::endl;
  return s;
}

/*******************************************************************/

std::ostream& operator<<(std::ostream &s, const ControlData& A)
{
  s << "ControlData\n";
  s << "status" <<"\t"<< A.status() << std::endl;
  s << "firstresidual" <<"\t"<< A.firstresidual() << std::endl;
  s << "residual" <<"\t"<< A.residual() << std::endl;
  s << "aimedresidual" <<"\t"<< A.aimedresidual() << std::endl;
  s << "previousresidual" << A.previousresidual() << std::endl;
  s << "iteration" <<"\t"<< A.iteration() << std::endl;
  s << "correction" <<"\t"<< A.correction() << std::endl;
  s << "previouscorrection" <<"\t"<< A.previouscorrection() << std::endl;
  return s;
}

/*******************************************************************/

std::ostream& operator<<(std::ostream &s, const UserData& A)
{
  s << "UserData\n";
  s << "tol" <<"\t"<< A.tol() << std::endl;
  s << "globaltol" <<"\t"<< A.globaltol() << std::endl;
  s << "breaktol" <<"\t"<< A.breaktol() << std::endl;
  s << "maxiter" <<"\t"<< A.maxiter() << std::endl;
  s << "printstep" <<"\t"<< A.printstep() << std::endl;
  return s;
}

/*******************************************************************/

std::ostream& operator<<(std::ostream &s, const CGInfo& A)
{
  s << "CGInfo\n";
  s << A.statistics()<<std::endl;
  s << A.control()<<std::endl;
  s << A.user()<<std::endl;
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
  double p = 1./GascoigneMath::max_int(1,CD.iteration());
  SD.rate() = pow(b,p);  
}

/*******************************************************************/

CGInfo::CGInfo(double f, double t, int p, int m, const std::string& txt)
{
  UD.text() = txt;
  
  UD.maxiter  () = m;
  UD.globaltol() = t;
  UD.printstep() = p;
  UD.tol()       = f;
  UD.breaktol()  = 1.e15;
 
  CD.reset();
  //SD.reset();
}

/*******************************************************************/

CGInfo::CGInfo(const std::string& txt)
{
  UD.text() = txt;
  double f = 1.e-6;
  double t=1.e-14;
  int p = 10; 
  int m = 100;
  
  UD.maxiter  () = m;
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
  if ( CD.residual()< GascoigneMath::max(CD.aimedresidual(),UD.globaltol()) )
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
      std::cout << UD.text() << " " << CD.iteration() << " Res: " << CD.residual();
      std::cout << " Cor: " << CD.correction()<< std::endl;
    }
  CD.iteration()++;
  if (CD.status()=="running")  return 0;

  return 1;
}
