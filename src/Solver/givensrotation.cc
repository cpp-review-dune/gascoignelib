#include "givensrotation.h"

using namespace Gascoigne;

/*--------------------------------------------------------------------------*/
 
GivensRotation::GivensRotation(int nn, double norm) : 
  n(0), H(nn,nn,0.), ci(nn-1,0.), si(nn-1,0.), gamma(nn,0.)
{ gamma[0] = norm;}

/*--------------------------------------------------------------------------*/
 
void GivensRotation::givens(double& a, double& b, int i) const
{
  double h = a;
  a =  ci[i]*h + si[i]*b;
  b = -si[i]*h + ci[i]*b;
}

/*--------------------------------------------------------------------------*/
 
double GivensRotation::orthogonalization(int dim) 
{
  assert(n==dim);
  int m = n;  n++;
  
  for (int i=0 ; i<m ; i++)
    {
      givens(H(i,m) ,H(i+1,m) ,i);
    }
  double beta = sqrt(H(m,m)*H(m,m) + H(n,m)*H(n,m));

  ci[m] = H(m,m) /beta;
  si[m] = H(n,m) /beta;
  gamma[n] = 0.;

  givens(H(m,m)     ,H(n,m)     ,m);
  givens(gamma[m],gamma[n] ,m);

  return gamma[n];
}

/*--------------------------------------------------------------------------*/
 
DoubleVector GivensRotation::getcoefficients()
{
  nvector<double> h(n);  
  for (int i=n-1; i>=0; i--)
    {
      double s = gamma[i];
      for (int j=i+1; j<n-1; j++) s -= H(i,j)*h[j];
      h[i] = s/H(i,i);
    }
  return h;
}

/*--------------------------------------------------------------------------*/
 
