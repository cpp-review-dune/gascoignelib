#include "baseq42dwithsecond.h"

#define NDOF   25
#define NDOF1d 5

/*------------------------------------------------------------------*/

namespace Gascoigne
{

BaseQ42dWithSecond::BaseQ42dWithSecond() : BaseQ42d()
{
  dxx.reservesize(NDOF);
  dyy.reservesize(NDOF);
  dxy.reservesize(NDOF);
}

/*------------------------------------------------------------------*/

void BaseQ42dWithSecond::point(const Vertex2d& s) const
{
  BaseQ42d::point(s);
  for(int i=0;i<NDOF;i++)
    {
      int ix = i%NDOF1d;
      int iy = i/NDOF1d;

      dxx[i]     = psi_xx(ix,s.x()) * psi   (iy,s.y());
      dxy[i]     = psi_x (ix,s.x()) * psi_x (iy,s.y());
      dyy[i]     = psi   (ix,s.x()) * psi_xx(iy,s.y());
    }
}

/*------------------------------------------------------------------*/

}

#undef NDOF
#undef NDOF1d
