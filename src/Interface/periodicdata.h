#ifndef  __periodicdata_h
#define  __periodicdata_h

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

    /// void operator()(Vertex2d& w, const Vertex2d& v, int col_w, int col_v)
    /// gets the coordinate v, color of boundarypart "col_v" and color 
    /// of boundarypart "col_w"and sets the coordinate w.
    ///
    //////////////////////////////////////////////

  class PeriodicData : public virtual Application
  {
    public:
      PeriodicData() : Application() {}

      virtual ~PeriodicData() {}

      virtual void operator()(DoubleVector& b, const Vertex2d& v, int col) const {
      /*-------------------------------------------------------
     | Falls auf einem der beiden Raender ein Offset addiert
     | werden soll, um einen Sprungterm in der Loesung zu
     | erhalten. Anwendung siehe DirichletData().
     | Nicht implementieren, falls echt periodische Loesung
     | gewuenscht ist.
     -------------------------------------------------------*/
      }

      virtual void operator()(DoubleVector& b, const Vertex3d& v, int col) const {
      /*-------------------------------------------------------
     | Falls auf einem der beiden Raender ein Offset addiert
     | werden soll, um einen Sprungterm in der Loesung zu
     | erhalten. Anwendung siehe DirichletData().
     | Nicht implementieren, falls echt periodische Loesung
     | gewuenscht ist.
     -------------------------------------------------------*/
      }
                  
      virtual std::set<int> preferred_colors()const {
        return std::set<int>();
      }
  };
}

#endif //__periodicdata_h

