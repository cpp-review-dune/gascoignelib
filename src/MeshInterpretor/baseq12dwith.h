#ifndef  __BaseQ12dWith_h
#define  __BaseQ12dWith_h

#include  "baseq12d.h"

/////////////////////////////////////////////
///
///@brief
///  computes in addition second derivatives

///
///
/////////////////////////////////////////////


class BaseQ12dWith : public BaseQ12d
{
private:

  nvector<double> dxy;

public:

  BaseQ12dWith();

  void   point(const Vertex2d& s);
  double phi_xy(int i) const {return dxy[i];}

};


#endif
