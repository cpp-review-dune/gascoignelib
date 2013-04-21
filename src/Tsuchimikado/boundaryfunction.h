/*----------------------------   boundaryfunction.h     ---------------------------*/
/*      $Id: boundaryfunction.h,v 1.1 2008/09/05 16:04:25 richter Exp $                 */
#ifndef __boundaryfunction_H
#define __boundaryfunction_H
/*----------------------------   boundaryfunction.h     ---------------------------*/

#include  <string>
#include  "vertex.h"

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
      virtual double operator()(const Gascoigne::Vertex<DIM>& c) const=0;
  
      virtual void grad(Gascoigne::Vertex<DIM>& dst, const Gascoigne::Vertex<DIM>& src) const;
      virtual void newton(Gascoigne::Vertex<DIM>&) const;
    };
}

/*---------------------------------------------------*/



/*----------------------------   boundaryfunction.h     ---------------------------*/
/* end of #ifndef __boundaryfunction_H */
#endif
/*----------------------------   boundaryfunction.h     ---------------------------*/
