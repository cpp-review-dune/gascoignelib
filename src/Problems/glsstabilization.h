#ifndef __GlsStabilization_h
#define __GlsStabilization_h

#include  "stabilization.h"

/*-------------------------------------------*/

class GlsStabilization : public Stabilization
{
 protected:

  double _delta, _tau;
  nvector<double> _sdelta;

  void NavierStokes(double h, double visc);

 public:

  GlsStabilization();

  double delta0, tau0, sdelta0;

  double  delta()      const { return _delta;}
  double& delta()            { return _delta;}
  double  delta(int i) const { return _sdelta[i];}
  double  tau()        const { return _tau;}
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
};

#endif
