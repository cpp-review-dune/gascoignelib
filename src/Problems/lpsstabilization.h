#ifndef __LpsStabilization_h
#define __LpsStabilization_h

#include  "stabilization.h"

/*-------------------------------------------*/

namespace Gascoigne
{
class LpsStabilization : public Stabilization
{
 protected:

  double _delta, _tau;
  nvector<double> _sdelta;

  void NavierStokes(double h, double visc);

 public:

  LpsStabilization();

  double delta0, sdelta0, tau0;

  double&  tau()      { return _tau;}
  double   tau()      const { return _tau;}
  double&  delta()      { return _delta;}
  double  delta()      const { return _delta;}
  double  delta(int i) const { return _sdelta[i];}
  void BasicInit(int n);
  void ReInit(double h, double visc, double u, double v)
    {
      norm(u,v);
      NavierStokes(h,visc);
    };
  void ReInit(double h, double visc, double u, double v, double w)
    {
      norm(u,v,w);
      NavierStokes(h,visc);
    };
  void ConvectionDiffusion(double visc);
  void ConvectionDiffusion(const nvector<double>& visc);
};
}

/*-------------------------------------------*/

#endif
