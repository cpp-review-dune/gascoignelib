#include  "runderkreis.h"

/*---------------------------------------------------*/

template<int DIM>
RunderKreis<DIM>::~RunderKreis() 
{
}

/*---------------------------------------------------*/

template<int DIM>
RunderKreis<DIM>::RunderKreis(double c, double sr) : 
  BoundaryFunction<DIM>(), 
  squareradius(sr) 
{ 
  center = c; 
}

/*---------------------------------------------------*/

template<int DIM>
RunderKreis<DIM>::RunderKreis(Vector c, double sr) : 
    BoundaryFunction<DIM>()
{
  center = c; 
  squareradius = sr;
}

/*---------------------------------------------------*/

template<int DIM>
Oval<DIM>::~Oval() {}

/*---------------------------------------------------*/

template<int DIM>
Oval<DIM>::Oval(const Vector& c, const Vector& sr) : center(c), radius(sr) {}

/*---------------------------------------------------*/

template RunderKreis<2>;
template RunderKreis<3>;
template Oval<2>;

