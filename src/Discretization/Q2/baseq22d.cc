#include "baseq22d.h"

#define NDOF   9
#define NDOF1d 3

/*------------------------------------------------------------------*/

namespace Gascoigne
{

BaseQ22d::BaseQ22d() 
{
  N.resize(NDOF);  
  DN.resize(NDOF);  

  //DDN.resize(NDOF);
  //for(int i=0;i<NDOF;i++) DDN[i].zero();

  a[0] = 1.;  b[0] = -3.;  c[0] = 2.;
  a[1] = 0.;  b[1] =  4.;  c[1] = -4.;
  a[2] = 0.;  b[2] = -1.;  c[2] = 2.;
}

/*------------------------------------------------------------------*/

void BaseQ22d::point(const Vertex2d& s) const
{
  for(int i=0;i<NDOF;i++) 
    { 
      int ix = i%NDOF1d;
      int iy = i/NDOF1d;

      N  [i]     = psi   (ix,s.x()) * psi   (iy,s.y());
      DN [i].x() = psi_x (ix,s.x()) * psi   (iy,s.y());
      DN [i].y() = psi   (ix,s.x()) * psi_x (iy,s.y());
//       DDN[i][0]  = psi_xx(ix,s.x()) * psi   (iy,s.y());
//       DDN[i][1]  = psi   (ix,s.x()) * psi_xx(iy,s.y());
//       DDN[i][2]  = psi_x (ix,s.x()) * psi_x (iy,s.y());
    }
}

/*------------------------------------------------------------------*/

}
#undef NDOF
#undef NDOF1d
