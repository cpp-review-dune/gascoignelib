#ifndef  __BaseQ23dWithSecond_h
#define  __BaseQ23dWithSecond_h


#include  "baseq23d.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////


class BaseQ23dWithSecond : public BaseQ23d
{
private:

  mutable DoubleVector dxx, dxy, dxz, dyy, dyz, dzz;

public:

  BaseQ23dWithSecond();

  void   point(const Vertex3d& s) const;

  double phi_xx(int i) const {return dxx[i];}
  double phi_xy(int i) const {return dxy[i];}
  double phi_xz(int i) const {return dxz[i];}
  double phi_yy(int i) const {return dyy[i];}
  double phi_yz(int i) const {return dyz[i];}
  double phi_zz(int i) const {return dzz[i];}
};
}


#endif
