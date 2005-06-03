#ifndef __IntegrationFormulaSummed_h
#define __IntegrationFormulaSummed_h

#include  "integrationformula.h"

/*------------------------------------------------------------*/

namespace Gascoigne
{

template<class INT>
class IntegrationFormulaSummed1d : public IntegrationFormula1d
{
 protected:

  INT I;

 public:

  IntegrationFormulaSummed1d(int n) : I(), IntegrationFormula1d()
    {
      int    N = static_cast<int>(pow(2.,n));

      IntegrationFormula1d::ReInit(N*I.n());

      int nn = static_cast<int>(pow(2.,n));
      double d2 = pow(0.5,n);

      for(int i=0;i<N;i++)
	{
	  int ix = i%nn;
	  double dx = 1.*ix;

	  for(int ii=0;ii<I.n();ii++)
	    {
	      int index = i*I.n()+ii;
	      w(index) = d2 * I.w(ii);

	      double x = I.c(ii).x();

	      x += dx;
	      x *= d2;

	      c(index).x() = x;
	    }
	}
    }
};

/*------------------------------------------------------------*/

template<class INT>
class IntegrationFormulaSummed2d : public IntegrationFormula2d
{
 protected:

  INT I;

 public:

  IntegrationFormulaSummed2d(int n=4) : I(), IntegrationFormula2d()
    {
      int    N = static_cast<int>(pow(4.,n));

      IntegrationFormula2d::ReInit(N*I.n());

      int nn = static_cast<int>(pow(2.,n));
      double d2 = pow(0.5,n);
      double d4 = pow(0.25,n);

      for(int i=0;i<N;i++)
	{
	  int ix = i%nn;
	  int iy = i/nn;

	  double dx = 1.*ix;
	  double dy = 1.*iy;

	  for(int ii=0;ii<I.n();ii++)
	    {
	      int index = i*I.n()+ii;
	      w(index) = d4 * I.w(ii);

	      double x = I.c(ii).x();
	      double y = I.c(ii).y();

	      x += dx;
	      y += dy;

	      x *= d2;
	      y *= d2;

	      c(index).x() = x;
	      c(index).y() = y;
	    }
	}
    }
};

}

#endif
