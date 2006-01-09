#ifndef  __BoundaryRightHandSide_h
#define  __BoundaryRightHandSide_h

#include  "vertex.h"
#include  <set>
#include  "nvector.h"
#include  <string>
#include  "gascoigne.h"
#include  "application.h"

/*-------------------------------------------------------*/

namespace Gascoigne
{
  
  //////////////////////////////////////////////
  ///
  ///@brief
  /// Interface class for Boundary Conditions of Neumann
  /// or Robin type

  /// void operator()(Vector& b, const Vertex2d& v, int col)
  /// gets the coordinate v and color of boundarypart "col" and 
  /// sets the values of b. b is a vector of length ncomp
  ///
  //////////////////////////////////////////////

  class BoundaryRightHandSide : public virtual Application
  {
    private:

    protected:

    public:
      BoundaryRightHandSide() : Application() {}
      ~BoundaryRightHandSide() {}

      virtual int GetNcomp() const=0;

      virtual void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v, const Vertex2d& n, int color) const {
        std::cerr << "\"BoundaryRightHandSide::operator()\" not written!" << std::endl;
        abort();
      }
      virtual void operator()(VectorIterator b, const TestFunction& N, const Vertex3d& v, const Vertex3d& n, int color) const {
        std::cerr << "\"BoundaryRightHandSide::operator()\" not written!" << std::endl;
        abort();
      }
  };

  typedef BoundaryRightHandSide BoundaryInitialCondition;

/*-------------------------------------------------------*/

}

#endif
