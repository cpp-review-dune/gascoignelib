#ifndef  __BaseQ22dWithSecond_h
#define  __BaseQ22dWithSecond_h

#include  "baseq22d.h"


/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////


class BaseQ22dWithSecond : public Gascoigne::BaseQ22d
{
private:
  
  mutable Gascoigne::nvector<double> dxx, dxy, dyy;

public:
  
  BaseQ22dWithSecond();

  void   point(const Gascoigne::Vertex2d& s) const;

  double phi_xx(int i) const {return dxx[i];}
  double phi_xy(int i) const {return dxy[i];}
  double phi_yy(int i) const {return dyy[i];}

};


#endif
