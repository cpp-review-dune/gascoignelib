#include "usefullfunctionsbd.h"
#include <math.h>

/**********************************************************/

namespace Gascoigne
{
double ParabelFunction(double x, double n0, double n1)
{
  return -4.*(x-n0)*(x-n1)/((n0-n1)*(n0-n1));
}

double GaussFunction(double x, double a, double stiff)
{
  return exp( -stiff*(x-a)*(x-a) );
}

double StiffnessFunction(double x, double a, double stiff)
{
  return 0.5 * (1.+tanh(stiff*(x-a)));
}

double PlugFlowFunction(double x, double a, double b, double stiff)
{
  double h = fabs(a-b);

  return tanh(stiff/h*(x-a)) * tanh(stiff/h*(-x+b));
}
}

/**********************************************************/

