#ifndef __DomainFunction_h
#define __DomainFunction_h

#include "application.h"
#include "vertex.h"

namespace Gascoigne
{

/**********************************************************/

class DomainFunction : public Application
{
  protected:

  public:
    DomainFunction() { }
    virtual ~DomainFunction() { }

    virtual int GetNcomp() const=0;
  
    virtual void F(DoubleVector& f, const Vertex2d &v) const
    {
      std::cerr << "DomainFunction::F not written\n";
      abort();
    }
    virtual void F(DoubleVector& f, const Vertex3d &v) const
    {
      std::cerr << "DomainFunction::F not written\n";
      abort(); 
    }
};

/**********************************************************/
}

#endif
