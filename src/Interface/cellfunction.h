#ifndef __CellFunction_h
#define __CellFunction_h

#include "application.h"
#include  "vertex.h"

namespace Gascoigne
{

/**********************************************************/

class CellFunction : public Application
{
  protected:

  public:
    CellFunction() : Application() { }
    ~CellFunction() { }

    virtual void F(LocalCellVector& b, double d) const
      {
	std::cout << "CellFunction::F not written\n";
	abort();
      }

    virtual void point(const Vertex2d& v) const {}
    virtual void point(const Vertex3d& v) const {}

    virtual int GetNcomp() const = 0;

};

/**********************************************************/
}

#endif
