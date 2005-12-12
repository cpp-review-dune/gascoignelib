#ifndef __CellRightHandSide_h
#define __CellRightHandSide_h

#include "application.h"
#include  "vertex.h"

namespace Gascoigne
{

/**********************************************************/

class CellRightHandSide : public Application
{
  protected:

  public:
    CellRightHandSide() : Application() { }
    ~CellRightHandSide() { }

    virtual void F(DoubleVector& b, double d) const
      {
	std::cout << "CellRightHandSide::F not written\n";
	abort();
      }

    virtual void point(const Vertex2d& v) const {}
    virtual void point(const Vertex3d& v) const {}

    virtual int GetNcomp() const = 0;

};

/**********************************************************/
}

#endif
