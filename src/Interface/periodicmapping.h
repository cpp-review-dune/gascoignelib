#ifndef  __periodicmapping_h
#define  __periodicmapping_h

#include  "vertex.h"
#include  <set>
#include  "nvector.h"
#include  <string>
#include  "application.h"


/*-----------------------------------------*/

namespace Gascoigne
{

    //////////////////////////////////////////////
    ///
    ///@brief
    /// Periodic Boundary Conditions

    /// void operator()(Vertex2d& w, const Vertex2d& v)
    /// gets the coordinate v of a vertex and sets the
    /// coordinate w of corresponding vertex on other
    /// boundary.
    ///
    //////////////////////////////////////////////

  class PeriodicMapping : public virtual Application
  {
    public:
      PeriodicMapping() : Application() {}

      virtual ~PeriodicMapping() {}

      virtual void transformCoords(Vertex2d& w, const Vertex2d& v) const = 0;
      virtual void transformCoords(Vertex3d& w, const Vertex3d& v) const = 0;

      virtual std::set<int> preferred_colors()const {
        return std::set<int>();
      }
  };
}

#endif //__periodicmapping_h

