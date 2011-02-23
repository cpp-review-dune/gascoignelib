#ifndef __baseq42d_h
#define __baseq42d_h

#include  <vector>
#include  <string>
#include  <utility>
#include  "vertex.h"
#include  "numfixarray.h"
#include  "base2d.h"

#define NDOF   25
#define NDOF1d 5

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  Basis on the reference element

///
///
/////////////////////////////////////////////

/**************************************************/

class BaseQ42d : public Base2d
{
 protected:

  fixarray<NDOF1d,double>     a,b,c,d,e;

  double psi   (int i, double x) const { return a[i] + b[i]*x + c[i]*x*x + d[i]*x*x*x + e[i]*x*x*x*x;}
  double psi_x (int i, double x) const { return b[i] + 2.*c[i]*x + 3.*d[i]*x*x + 4.*e[i]*x*x*x;}
  double psi_xx(int i, double x) const { return 2.*c[i] + 6.*d[i]*x + 12.*e[i]*x*x; }

 public:

  BaseQ42d();

  int  n()             const { return NDOF;}
  void point(const Vertex2d& s) const;

  double phi   (int i) const { return N  [i];}
  double phi_x (int i) const { return DN [i].x();}
  double phi_y (int i) const { return DN [i].y();}
  double phi_xx(int i) const {
    std::cerr << "\"BaseQ42d::phi_xx\" not written!" << std::endl;
    abort();
  }
  double phi_yy(int i) const {
    std::cerr << "\"BaseQ42d::phi_yy\" not written!" << std::endl;
    abort();
  }
  double phi_xy(int i) const {
    std::cerr << "\"BaseQ42d::phi_xy\" not written!" << std::endl;
    abort();
  }

  const Vertex2d&  phi_grad (int i) const {return DN[i];}
};
}

#undef NDOF
#undef NDOF1d

#endif
