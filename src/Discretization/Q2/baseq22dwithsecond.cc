#include  "baseq22dwithsecond.h"


using namespace std;
using namespace Gascoigne;

#define NDOF   9
#define NDOF1d 3

/*-----------------------------------------*/

BaseQ22dWithSecond::BaseQ22dWithSecond() : BaseQ22d()
{
  dxx.reservesize(NDOF);
  dxy.reservesize(NDOF);
  dyy.reservesize(NDOF);
}

/*------------------------------------------------------------------*/

void BaseQ22dWithSecond::point(const Vertex2d& s) const 
{
  BaseQ22d::point(s);
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

#undef NDOF
#undef NDOF1d
