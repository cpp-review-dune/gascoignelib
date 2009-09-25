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
                  
      virtual void transformCoords(Vertex2d& w, const Vertex2d& v, int col_w, int col_v) const {
      /*-------------------------------------------------------
     | Affin-lineare Transformation, die die beiden
     | periodischen Raender aufeinander abbildet.
     | v = gegebener Knoten, w = gesuchter Knoten auf dem
     | anderen Rand.
     | Die hier angegebene Identitaet funktioniert nur, wenn
     | es sich um eine Translation im rechten Winkel zum Rand
     | handelt und der Rand ein Geradenstueck ist.
     | Sonst lokal wie gewuenscht ueberschreiben.
     -------------------------------------------------------*/
        w.x() = v.x();
        w.y() = v.y();
      }
      
      virtual void transformCoords(Vertex3d& w, const Vertex3d& v, int col_w, int col_v) const {
      /*-------------------------------------------------------
     | Affin-lineare Transformation, die die beiden
     | periodischen Raender aufeinander abbildet.
     | v = gegebener Knoten, w = gesuchter Knoten auf dem
     | anderen Rand.
     | Die hier angegebene Identitaet funktioniert nur, wenn
     | es sich um eine Translation im rechten Winkel zum Rand
     | handelt und der Rand ein Geradenstueck ist.
     | Sonst lokal wie gewuenscht ueberschreiben.
     -------------------------------------------------------*/
        w.x() = v.x();
        w.y() = w.y();
        w.z() = w.z();
      }
      
      virtual std::set<int> preferred_colors()const {
        return std::set<int>();
      }
  };
}

#endif //__periodicdata_h

