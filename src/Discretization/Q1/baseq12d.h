#ifndef __baseq12d_h
#define __baseq12d_h

#include  <vector>
#include  <string>
#include  <utility>
#include  "vertex.h"
#include  "numfixarray.h"
#include  "base2d.h"

/////////////////////////////////////////////
///
///@brief
///  Basis on the reference element

///
///
/////////////////////////////////////////////

#define NDOF   4
#define NDOF1d 2

/**************************************************/

namespace Gascoigne
{
class BaseQ12d : public Base2d
{
 protected:

  fixarray<2,double>          a,b;
  mutable DoubleVector dxy;

  void BasicInit();

  double psi_x(int i, double x) const { return b[i]; }

 public:
  
  BaseQ12d();

  double psi  (int i, double x) const { return a[i] + b[i]*x;}
  int    n()                    const { return NDOF;}
  void   point(const Vertex2d& s) const;

  double phi   (int i) const {return N  [i];}
  double phi_x (int i) const {return DN [i].x();}
  double phi_y (int i) const {return DN [i].y();}
  double phi_xx(int i) const {return 0.;}
  double phi_yy(int i) const {return 0.;}
/*   double phi_xy(int i) const {return 0.;} */
  double phi_xy(int i) const {assert(0); return dxy[i];}
  const Vertex2d &  phi_grad (int i) const {return DN [i];}
};
}

#undef NDOF
#undef NDOF1d

#endif
