#ifndef __vertex_h
#define __vertex_h

#include "numfixarray.h"

/*------------------------------------------------------------*/
namespace Gascoigne
{

template<int N>
class Vertex : public numfixarray<N,double>
{
public:

  Vertex<N>() : numfixarray<N,double>() {}
  Vertex<N>(const Vertex& c) : numfixarray<N,double>(c) {}
  Vertex<N>(const double& x0) : numfixarray<N,double>(x0) {}
  Vertex<N>(const double& x0, const double& y0) : numfixarray<N,double>() 
    {x()=x0; y()=y0;}
  Vertex<N>(const double& x0, const double& y0, const double& z0) : numfixarray<N,double>() 
    {x()=x0; y()=y0; z()=z0;}
 
  Vertex<N>& operator=(const Vertex<N>& c) 
    {
      numfixarray<N,double>::operator=(c);
      return *this;
    }
  Vertex<N>& operator=(double d) 
    {
      numfixarray<N,double>::operator=(d);
      return *this;
    }

  const double& x() const  { return (*this)[0]; }
  const double& y() const  { return (*this)[1]; }
  const double& z() const  { return (*this)[2]; }
  double&       x()        { return (*this)[0]; }
  double&       y()        { return (*this)[1]; }
  double&       z()        { return (*this)[2]; }
};

/*------------------------------------------------------------*/

typedef  Vertex<1>   Vertex1d; 
typedef  Vertex<2>   Vertex2d; 
typedef  Vertex<3>   Vertex3d; 
}

#endif
