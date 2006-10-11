#ifndef __DomainRightHandSide_h
#define __DomainRightHandSide_h

#include  "application.h"
#include  "vertex.h"

namespace Gascoigne
{

/*-------------------------------------------------------*/

  class DomainRightHandSide : public virtual Application
  {
    private:

    protected:

    public:
      DomainRightHandSide() { }
      ~DomainRightHandSide() { }

      virtual int GetNcomp() const=0;

      virtual double operator()(int c, const Vertex2d& v) const {
        std::cerr << "\"DomainRightHandSide::operator()\" not written" << std::endl;
        abort();
      }
      virtual double operator()(int c, const Vertex3d& v) const {
        std::cerr << "\"DomainRightHandSide::operator()\" not written" << std::endl;
        abort();
      }

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

      virtual void SetCellSize(double h) const { }
  };
  
  typedef DomainRightHandSide DomainInitialCondition;

/*-------------------------------------------------------*/

}

#endif
