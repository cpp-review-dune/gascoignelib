#ifndef __stabilization_h
#define __stabilization_h

#include "nvector.h"

/*-------------------------------------------*/

class Stabilization
{
 protected:

  double _alpha, _h, dt, _dtfactor, _norm;

  void norm(double u, double v)
    { 
      _norm = sqrt(u*u+v*v)+1e-6; 
    }
  void norm(double u, double v, double w) 
    { 
      _norm = sqrt(u*u+v*v+w*w)+1e-6; 
    }

 public:

  Stabilization() :  xeta0(6.), alpha0(0.2), dt(0.), _dtfactor(1.), _norm(0.) {}
  ~Stabilization() {}

  double alpha0, xeta0;

  double& DeltaT()           { return dt;}
  double  DeltaT()     const { return dt;}
  double  alpha()      const { return _alpha;}
  double& alpha()            { return _alpha;}
  double  h()          const { return _h;}
  double& dtfactor()         { return _dtfactor; }
  void    SetConvection(double u, double v)           { norm(u,v);}    
  void    SetConvection(double u, double v, double w) { norm(u,v,w);}    

  void ReInit(double h, double visc)
    {
      _h = h;
      double denom = xeta0 * visc / (h*h);
      if(dt>0.)
	{
	  denom += _dtfactor/dt;
	}
      _alpha = alpha0 / denom;
    }
};

#endif
