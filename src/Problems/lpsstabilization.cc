#include  "lpsstabilization.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{
LpsStabilization::LpsStabilization() : Stabilization(), _sdelta(0)
{
  _delta = delta0 = sdelta0 = 0.;
}

/*--------------------------------------------------------------*/

void LpsStabilization::BasicInit(int n)
{
  _sdelta.resize(n,0.);
}

/*--------------------------------------------------------------*/

void LpsStabilization::NavierStokes(double h, double visc)
{
  _h = h;

  double val  = xeta0 * visc/(h*h);
  double valc = xeta0 * visc/(h*h) + _norm/h;
  assert(dt>0);
  if(dt>0.)
    {
      val  += _dtfactor/dt;
      valc += _dtfactor/dt;
    }
  _alpha = alpha0 / val;
  _delta = delta0 / valc;
}

/*--------------------------------------------------------------*/

void LpsStabilization::ConvectionDiffusion(double visc)
{
  double val = xeta0 * visc/(_h*_h) + _norm/_h;
  if(dt>0.)
    {
      val  += _dtfactor/dt;
    }
  _sdelta[0] = sdelta0 / val; 
}

/*--------------------------------------------------------------*/

void LpsStabilization::ConvectionDiffusion(const nvector<double>& visc)
{
  double val = _norm/_h;
  if(dt>0.)
    {
      val += _dtfactor/dt;
    }
  for (int i=0; i<_sdelta.size(); i++)
    {
      _sdelta[i] = sdelta0 / (val+xeta0 * visc[i]/(_h*_h) ); 
    }
}
}
