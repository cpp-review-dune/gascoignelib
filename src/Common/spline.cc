#include "spline.h"

/*******************************************************/

CubicSpline::CubicSpline(const Vector& xx, const Vector& yy) :
  xa(xx), y(yy), m(xx.size())
{
  n = xx.size();

  yp0 = 0.; // df(x_0)
  yp1 = 0.; // df(x_n)

  compute_moments();
}

/*******************************************************/

double CubicSpline::lambda(int i) const
{
  if (i==0)     return 1.;
  return h(i+1)/(h(i)+h(i+1));
}

/*******************************************************/

double CubicSpline::mu(int i) const
{
  if (i==0) return 1.;
  return 1.-lambda(i);
}

/*******************************************************/

double CubicSpline::dcoeff(int i) const
{
  if (i==0)   return  6*(df(0)  -yp0)/h(1);
  if (i==n-1) return -6*(df(n-1)-yp1)/h(n-1);
  return 6.*ddf(i);
}

/*******************************************************/

void CubicSpline::compute_moments()
{
  nvector<double> u(n,0.), q(n,0.);

  q[0] = -lambda(0)/2.;
  u[0] =  dcoeff(0) /2.;

  for (int k=1; k<n; k++)
    {
      double p = mu(k)*q[k-1]+2.;
      q[k] = -lambda(k)/p;
      u[k] = (dcoeff(k)-mu(k)*u[k-1])/p;
    }
  m[n-1] = u[n-1];
  for (int i=n-2; i>=0; i--)
    {
      m[i] = q[i]*m[i+1]+u[i];
    }
}

/*******************************************************/

double CubicSpline::operator()(double x) const
{
  for (int i=0; i<n-1; i++)
    {
      if ((x>=xa[i]) && (x<=xa[i+1]))
	{
	  double alpha = y[i];
	  double beta  = df(i) - (2.*m[i]+m[i+1]) * h(i+1)/6.;
	  double gamma = m[i]*0.5;
	  double delta = (m[i+1]-m[i])/(6.*h(i+1));
	  double d = x-xa[i];

	  return alpha + beta*d + gamma*d*d + delta*d*d*d;
	}
    }
  std::cout << "error: spline out of range " << x << std::endl;
  exit(1);
  return 0.;
}

