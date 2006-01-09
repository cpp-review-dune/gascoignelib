#ifndef  __BoundaryFunctional_h
#define  __BoundaryFunctional_h

#include  <set>
#include  "functional.h"
#include  "vertex.h"

/*-----------------------------------------*/


namespace Gascoigne
{
  class BoundaryFunctional : public virtual Functional
  {
    private:

    protected:

    public:
      BoundaryFunctional() : Functional() {}
      virtual ~BoundaryFunctional() {};

      virtual double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const {
        std::cerr << "\"BoundaryFunctional::J\" for 2d not written!" << std::endl;
        abort();
      }
      virtual double J(const FemFunction& U, const Vertex3d& v, const Vertex3d& n, int color) const {
        std::cerr << "\"BoundaryFunctional::J\" for 3d not written!" << std::endl;
        abort();
      }
  };
}

#endif
