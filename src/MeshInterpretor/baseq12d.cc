#include "baseq12d.h"

#define NDOF   4
#define NDOF1d 2

/*------------------------------------------------------------------*/

namespace Gascoigne
{
BaseQ12d::BaseQ12d()  
{
  BasicInit();
}

/*------------------------------------------------------------------*/

void BaseQ12d::BasicInit()
{
  N.resize(NDOF);  
  DN.resize(NDOF);  
//   dxy.reservesize(NDOF);

  a[0] = 1.;  b[0] = -1.;
  a[1] = 0.;  b[1] =  1.;
}

/*------------------------------------------------------------------*/

void BaseQ12d::point(const Vertex2d& s) const
{
  for(int i=0;i<NDOF;i++) 
    { 
      int ix = i%NDOF1d;
      int iy = i/NDOF1d;

      N  [i]     = psi   (ix,s.x()) * psi   (iy,s.y());
      DN [i].x() = psi_x (ix,s.x()) * psi   (iy,s.y());
      DN [i].y() = psi   (ix,s.x()) * psi_x (iy,s.y());
//       dxy[i]     = psi_x(ix,s.x()) * psi_x(iy,s.y());
    }
}
}

/*------------------------------------------------------------------*/

#undef NDOF
#undef NDOF1d
