#ifndef __spline_h
#define __spline_h

#include "nvector.h"

/*******************************************************/

class CubicSpline
{
  protected:
  
  typedef nvector<double>  Vector;

  Vector xa, y, m;
  int n;
  double yp0, yp1;

  void compute_moments();
  double h  (int i) const { return   xa[i] - xa[i-1]; }
  double df (int i) const { return ( y[i+1] - y[i] ) / h(i+1); }
  double ddf(int i) const { return ( df(i) - df(i-1) ) / (h(i+1)+h(i)); }
  double lambda(int i) const;
  double mu    (int i) const;
  double dcoeff(int i) const;
  
  public:

  CubicSpline(const Vector&, const Vector&); 
  double operator()(double x) const;
};

/*******************************************************/

#endif
