#ifndef __baseq43dwithsecond_h
#define __baseq43dwithsecond_h

#include  "baseq43d.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////

class BaseQ43dWithSecond : public BaseQ43d
{

 protected:

  mutable DoubleVector dxx,dyy,dzz,dxy,dxz,dyz;

 public:

  BaseQ43dWithSecond();

  void point(const Vertex3d& s) const;

  double phi_xx(int i) const {return dxx[i];}
  double phi_yy(int i) const {return dyy[i];}
  double phi_zz(int i) const {return dzz[i];}
  double phi_xy(int i) const {return dxy[i];}
  double phi_xz(int i) const {return dxz[i];}
  double phi_yz(int i) const {return dyz[i];}
};
}

#endif
