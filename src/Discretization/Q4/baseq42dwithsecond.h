#ifndef __baseq42dwithsecond_h
#define __baseq42dwithsecond_h

#include  "baseq42d.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////

/**************************************************/

class BaseQ42dWithSecond : public BaseQ42d
{
 protected:

  mutable DoubleVector dxx, dxy, dyy;

 public:

  BaseQ42dWithSecond();

  void point(const Vertex2d& s) const;

  double phi_xx(int i) const { return dxx[i];}
  double phi_yy(int i) const { return dyy[i];}
  double phi_xy(int i) const { return dxy[i];}
};
}

#endif
