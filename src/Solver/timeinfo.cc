#include  "timeinfo.h"

#include  <iostream>
#include  <stdio.h>
#include  <cassert>
#include  <math.h>

using namespace std;

/*-----------------------------------------*/

TimeInfo::TimeInfo()
{
  BasicInit();
}

/*-----------------------------------------*/

void TimeInfo::BasicInit()
{
  // fractional theta parameter
  //
  double gamma = 1.-sqrt(0.5);

  _ftscale[0] = gamma;
  _ftscale[1] = 1.-2.*gamma;
  _ftscale[2] = gamma;

  double alpha = (1.-2*gamma)/(1.-gamma);
  double beta = 1.-alpha;

  _fttheta[0] = alpha;
  _fttheta[1] = beta;
  _fttheta[2] = alpha;

  _theta = 1.;
  _scheme = "Euler";

  ReInit();
}

/*-----------------------------------------*/

void TimeInfo::ReInit()
{
  _time = 0.;
  _tbegin = _tend = _deltat = 0.;
  _iter = 0;
}

/*-----------------------------------------*/

void TimeInfo::ReInit(double det)
{
  _time   = _tbegin;
  _deltat = det;
  _iter = 0;
}

/*-----------------------------------------*/

int TimeInfo::ftstep() const
{
  if (_iter<=1) return 0;
  return (_iter-1)%3;
}

/*-----------------------------------------*/

void TimeInfo::ReInit(double tb, double te, double det, const string& sch, int ne, double tt)
{
  _tbegin = tb;
  _time      = _tbegin;
  _tend   = te;
  _deltat = det;
  _neuler = ne;
  _Theta  = tt;

  _scheme = sch;
  if      (_scheme=="CN"   )           _theta = 0.5;
  else if (_scheme=="Euler")           _theta = 1.;
  else if (_scheme=="FractionalTheta") _theta = 0.;
  else if (_scheme=="Theta")           _theta = _Theta;
  else assert(0);

  if (_scheme=="FractionalTheta") _deltat *= 3.;
}

/*-----------------------------------------*/

double TimeInfo::theta() const 
{ 
  if (_iter<=_neuler) return 1.;
  if (_scheme=="FractionalTheta") return _fttheta[ftstep()];
  return _theta;
}

/*-----------------------------------------*/

double TimeInfo::dt() const 
{
  double h = 1.;
  if ((_iter>_neuler)&&(_scheme=="FractionalTheta")) h = _ftscale[ftstep()];
  
  assert((_deltat==0) || (_deltat>1.e-10));

  return _deltat * h;
}

/*-----------------------------------------*/

void TimeInfo::ReInitBackward(int niter, double endtime)
{
  _iter = niter;
  _time = endtime;
}

/*-----------------------------------------*/

void TimeInfo::iteration_backward(int i)
{
  assert(i<=_iter);

  _iter = i;
  _time -= dt();

  string actualscheme = _scheme;
  if (_iter<=_neuler) actualscheme = "Euler";

  cout << "\n============= " << _iter << " ========== " << actualscheme;
  cout << " [t,dt] "<< _time << " " << dt() << "\t";
  cout << endl;
}

/*-----------------------------------------*/

void TimeInfo::iteration(int i)
{
  assert(i>=_iter);

  _iter = i;
  _time += dt();

  string actualscheme = _scheme;
  if (_iter<=_neuler) actualscheme = "Euler";

  cout << "\n============= " << _iter << " ========== " << actualscheme;
  cout << " [t,dt] "<< _time << " " << dt() << "\t";
  cout << endl;
}

/*-----------------------------------------*/

