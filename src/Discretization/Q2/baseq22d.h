#ifndef __baseq22d_h
#define __baseq22d_h

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

#define NDOF   9
#define NDOF1d 3

namespace Gascoigne
{
/**************************************************/

class BaseQ22d : public Base2d
{
 protected:

  fixarray<NDOF1d,double>     a,b,c;
  
  double psi   (int i, double x) const { return a[i] + b[i]*x + c[i]*x*x;}
  double psi_x (int i, double x) const { return b[i] + 2.*c[i]*x;}
  double psi_xx(int i, double x) const { return 2.*c[i]; }

 public:
  
  BaseQ22d();

  int  n()             const { return NDOF;}
  void point(const Vertex2d& s) const;

  double phi   (int i) const { return N  [i];}
  double phi_x (int i) const { return DN [i].x();}
  double phi_y (int i) const { return DN [i].y();}
  double phi_xx(int i) const { assert(0);}
  double phi_yy(int i) const { assert(0);}
  double phi_xy(int i) const { assert(0);}

  const Vertex2d&  phi_grad (int i) const {return DN[i];}
};
}

#undef NDOF
#undef NDOF1d

#endif
