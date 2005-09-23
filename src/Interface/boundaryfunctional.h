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

      virtual std::set<int> GetColors() const=0;
      virtual double J(const FemFunction& U, const Vertex3d& v) const {
        std::cerr << "\"BoundaryFunctional::J\" for 3d not written!" << std::endl;
        abort();
      }
      virtual double J(const FemFunction& U, const Vertex2d& v) const {
        std::cerr << "\"BoundaryFunctional::J\" for 2d not written!" << std::endl;
        abort();
      }
      virtual void J(DoubleVector& b, const FemFunction& U, const TestFunction& N) const {
        std::cerr << "\"BoundaryFunctional::J\" not written!" << std::endl;
        abort();
      }
  };
}

#endif
