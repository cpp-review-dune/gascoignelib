#ifndef  __BaseQ23dWithSecond_h
#define  __BaseQ23dWithSecond_h


#include  "baseq23d.h"


/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////


class BaseQ23dWithSecond : public Gascoigne::BaseQ23d
{
private:

  mutable Gascoigne::nvector<double> dxx, dxy, dxz, dyy, dyz, dzz;

public:

  BaseQ23dWithSecond();

  void   point(const Gascoigne::Vertex3d& s) const;

  double phi_xx(int i) const {return dxx[i];}
  double phi_xy(int i) const {return dxy[i];}
  double phi_xz(int i) const {return dxz[i];}
  double phi_yy(int i) const {return dyy[i];}
  double phi_yz(int i) const {return dyz[i];}
  double phi_zz(int i) const {return dzz[i];}

};


#endif
