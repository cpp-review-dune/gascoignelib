#include  "glsstabilization.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{
GlsStabilization::GlsStabilization() : Stabilization()
{
  _delta = _tau = 0.;
  delta0 = sdelta0 = tau0 = 0.;
  _sdelta.resize(1,0.);
}

/*--------------------------------------------------------------*/

void GlsStabilization::NavierStokes(double h, double visc)
{
  _h = h;

  double val  = xeta0 * visc/(h*h) + _norm/h;
  double valc = xeta0 * visc/(h*h) + _norm/h;
  if(dt>0.)
    {
      valc += _dtfactor/dt;
    }
  _alpha = alpha0 / val;
  _delta = delta0 / valc;
  _tau   = tau0   * _norm * _norm *_delta;
}

/*--------------------------------------------------------------*/

void GlsStabilization::ConvectionDiffusion(double visc)
{
  double val = xeta0 * visc/(_h*_h) + _norm/_h;
  if(dt>0.)
    {
      val  += _dtfactor/dt;
    }
  _sdelta[0] = sdelta0 / val;
}
}
