#ifndef  __BaseQ13dWith_h
#define  __BaseQ13dWith_h

#include  "baseq13d.h"

/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////


class BaseQ13dWith : public BaseQ13d
{
private:

  nvector<double> dxy, dxz, dyz;

public:

  BaseQ13dWith();

  void   point(const Vertex3d& s);
  double phi_xy(int i) const {return dxy[i];}
  double phi_xz(int i) const {return dxz[i];}
  double phi_yy(int i) const {return dyz[i];}

};


#endif
