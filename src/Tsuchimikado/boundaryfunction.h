/*----------------------------   boundaryfunction.h     ---------------------------*/
/*      $Id: boundaryfunction.h,v 1.1 2008/09/05 16:04:25 richter Exp $                 */
#ifndef __boundaryfunction_H
#define __boundaryfunction_H
/*----------------------------   boundaryfunction.h     ---------------------------*/

#include  <string>
#include  "meshvertex.h"

/*---------------------------------------------------*/

namespace Tsuchimikado
{
  template<int DIM>
    class BoundaryFunction
    {
    protected:
  
    public :
  
      BoundaryFunction() {}
      virtual ~BoundaryFunction() {}

      virtual std::string GetName() const=0;
      virtual double operator()(const MeshVertex<DIM>& c) const=0;
  
      virtual void grad(MeshVertex<DIM>& dst, const MeshVertex<DIM>& src) const;
      virtual void newton(MeshVertex<DIM>&) const;
    };
}

/*---------------------------------------------------*/



/*----------------------------   boundaryfunction.h     ---------------------------*/
/* end of #ifndef __boundaryfunction_H */
#endif
/*----------------------------   boundaryfunction.h     ---------------------------*/
