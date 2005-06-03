#ifndef __CellFunction_h
#define __CellFunction_h

#include "application.h"
#include  "vertex.h"

namespace Gascoigne
{

/**********************************************************/

class CellFunction : public Gascoigne::Application
{
  protected:

  public:
    CellFunction() : Gascoigne::Application() { }
    ~CellFunction() { }

    virtual void F(Gascoigne::LocalCellVector& b, double d) const
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
