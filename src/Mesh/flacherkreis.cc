#include  "runderkreis.h"

/*---------------------------------------------------*/

FlacherKreis::~FlacherKreis() {}

/*---------------------------------------------------*/

FlacherKreis::FlacherKreis(double c, double sr) : 
    BoundaryFunction<3>(), 
    squareradius(sr) 
{ 
  center = c; 
}

/*---------------------------------------------------*/

FlacherKreis::FlacherKreis(Vector& c, double sr) : 
  BoundaryFunction<3>(), 
  center(c), 
  squareradius(sr) 
{}

/*---------------------------------------------------*/

double FlacherKreis::operator()(const Vector& c) const 
{
  double r = - squareradius;
  for (int i=0; i<2; i++)
    {
      double dx = c[i]-center[i];
      r += dx * dx;
    }
  return r;
}
