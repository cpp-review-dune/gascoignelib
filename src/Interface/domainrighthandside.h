#ifndef __DomainRightHandSide_h
#define __DomainRightHandSide_h

#include "righthandsidedata.h"

namespace Gascoigne
{

/**********************************************************/

  class DomainRightHandSide : public virtual RightHandSideData
  {
    private:

    protected:

    public:
      DomainRightHandSide() : RightHandSideData() { }
      ~DomainRightHandSide() { }

      virtual double operator()(int c, const Vertex2d& v) const {assert(0); return 0;}
      virtual double operator()(int c, const Vertex3d& v) const {assert(0); return 0;}

      virtual void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
        {
          for(int c=0;c<GetNcomp();c++)
            {
              b[c] += N.m()* (*this)(c,v);
            }
        }
      virtual void operator()(VectorIterator b, const TestFunction& N, const Vertex3d& v) const 
        {
          for(int c=0;c<GetNcomp();c++)
            {
              b[c] += N.m()* (*this)(c,v);
            }
        }
  };

/**********************************************************/

}

#endif
