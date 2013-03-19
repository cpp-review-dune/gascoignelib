/*----------------------------   meshvertex.h     ---------------------------*/
/*      $Id: meshvertex.h,v 1.1.1.1 2008/09/03 21:06:14 richter Exp $                 */
#ifndef __meshvertex_H
#define __meshvertex_H
/*----------------------------   meshvertex.h     ---------------------------*/


#include "mesharray.h"

/*------------------------------------------------------------*/
namespace Tsuchimikado
{

  /**
   * simple structure for a vertex used in the triacontainer
   **/
  template<int N>
    class MeshVertex : public mesharray<N,double>
  {
  public:

    MeshVertex<N>() : mesharray<N,double>() {}
    MeshVertex<N>(const MeshVertex& c) : mesharray<N,double>(c) {}
    MeshVertex<N>(const double& x0) : mesharray<N,double>(x0) {}
    MeshVertex<N>(const double& x0, const double& y0) : mesharray<N,double>() 
    {x()=x0; y()=y0;}
    MeshVertex<N>(const double& x0, const double& y0, const double& z0) : mesharray<N,double>() 
    {x()=x0; y()=y0; z()=z0;}
 
    MeshVertex<N>& operator=(const MeshVertex<N>& c) 
      {
	mesharray<N,double>::operator=(c);
	return *this;
      }
    MeshVertex<N>& operator=(double d) 
      {
	mesharray<N,double>::operator=(d);
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

  typedef  MeshVertex<1>   MeshVertex1d; 
  typedef  MeshVertex<2>   MeshVertex2d; 
  typedef  MeshVertex<3>   MeshVertex3d; 
}



/*----------------------------   meshvertex.h     --------------------------\-*/
/* end of #ifndef __meshvertex_H */
#endif
/*----------------------------   meshvertex.h     ---------------------------*/
