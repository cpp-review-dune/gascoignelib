#ifndef __IntegrationFormulaSummed_h
#define __IntegrationFormulaSummed_h

#include  "integrationformula.h"

/*------------------------------------------------------------*/

template<class INT>
class IntegrationFormulaSummed1d : public IntegrationFormula1d
{
 protected:

  INT I;

 public:

  IntegrationFormulaSummed1d(int n) : I(), IntegrationFormula1d()
    {
      int    N = (int) pow(2,n);

      IntegrationFormula1d::init(N*I.n());

      int nn = (int) pow(2,n);
      double d2 = pow(0.5,n);

      for(int i=0;i<N;i++)
	{
	  int ix = i%nn;
	  double dx = 1.*ix;

	  for(int ii=0;ii<I.n();ii++)
	    {
	      int index = i*I.n()+ii;
	      iw[index] = d2 * I.w(ii);

	      double x = I.xi(ii).x();

	      x += dx;
	      x *= d2;

	      ic[index].x() = x;
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

  IntegrationFormulaSummed2d(int n) : I(), IntegrationFormula2d()
    {
      int    N = (int) pow(4,n);

      IntegrationFormula2d::init(N*I.n());

      int nn = (int) pow(2,n);
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
	      iw[index] = d4 * I.w(ii);

	      double x = I.xi(ii).x();
	      double y = I.xi(ii).y();

	      x += dx;
	      y += dy;

	      x *= d2;
	      y *= d2;

	      ic[index].x() = x;
	      ic[index].y() = y;
	    }
	}
    }
};

#endif
