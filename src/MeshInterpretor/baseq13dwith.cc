#include  "baseq13dwith.h"

using namespace std;

#define NDOF   8
#define NDOF1d 2

/*-----------------------------------------*/

BaseQ13dWith::BaseQ13dWith() : BaseQ13d()
{
  dxy.reservesize(NDOF);
  dxz.reservesize(NDOF);
  dyz.reservesize(NDOF);
}

/*------------------------------------------------------------------*/

void BaseQ13dWith::point(const Vertex3d& s)
{
  BaseQ13d::point(s);
  for(int i=0;i<NDOF;i++) 
    { 
      int ix = i%2;
      int iy = (i%4)/2;
      int iz = i/4;

      dxy[i]     = psi_x(ix,s.x()) * psi_x(iy,s.y()) * psi  (iz,s.z());
      dxz[i]     = psi_x(ix,s.x()) * psi  (iy,s.y()) * psi_x(iz,s.z());
      dyz[i]     = psi  (ix,s.x()) * psi_x(iy,s.y()) * psi_x(iz,s.z());
    }
}

/*------------------------------------------------------------------*/

#undef NDOF
#undef NDOF1d
