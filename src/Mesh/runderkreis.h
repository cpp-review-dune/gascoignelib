#ifndef  __runderkreis_h
#define  __runderkreis_h

#include  "boundaryfunction.h"

/*---------------------------------------------------*/

template<int DIM>
class RunderKreis : public BoundaryFunction<DIM>
{
  typedef typename BoundaryFunction<DIM>::Vector Vector;

  double squareradius;
  Vector center;

public :

  ~RunderKreis();
  RunderKreis(double c, double sr);
  RunderKreis(Vector c, double sr);

  double operator()(const Vector& c) const 
    {
      double r = - squareradius;
      for (int i=0; i<DIM; i++)
	{
	  double dx = c[i]-center[i];
	  r += dx * dx;
	}
      return r;
    }
};

/*---------------------------------------------------*/

template<int DIM>
class Oval : public BoundaryFunction<DIM>
{
  typedef typename BoundaryFunction<DIM>::Vector Vector;

  Vector center, radius;

public :

  ~Oval();
  Oval(const Vector& c, const Vector& sr);

  double operator()(const Vector& c) const 
    {
      double r = -1.;
      for (int i=0; i<DIM; i++)
	{
	  double dx = (c[i]-center[i])/radius[i];
	  r += dx * dx;
	}
      return r;
    }
};

/*---------------------------------------------------*/

class FlacherKreis : public BoundaryFunction<3>
{
  typedef numfixarray<3,double>  Vector;

  double  squareradius;
  Vector  center;

public :

   ~FlacherKreis();
  FlacherKreis(double c, double sr);
  FlacherKreis(Vector& c, double sr);

  double operator()(const Vector& c) const;
};

/*---------------------------------------------------*/

#endif
