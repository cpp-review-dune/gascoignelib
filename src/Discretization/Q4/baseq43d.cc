#include "baseq43d.h"

#define NDOF   125
#define NDOF1d 5
#define NDOF2d 25

namespace Gascoigne
{

/*------------------------------------------------------------------*/

BaseQ43d::BaseQ43d()
{
  N  .resize(NDOF);
  DN .resize(NDOF);

  a[0] = 1.; b[0] = -25./3; c[0] = 70./3;   d[0] = -80./3;  e[0] = 32./3;
  a[1] = 0;  b[1] = 16.;    c[1] = -208./3; d[1] = 96.;     e[1] = -128./3;
  a[2] = 0;  b[2] = -12.;   c[2] = 76.;     d[2] = -128.;   e[2] = 64.;
  a[3] = 0;  b[3] = 16./3;  c[3] = -112./3; d[3] = 224./3.; e[3] = -128./3.;
  a[4] = 0;  b[4] = -1;     c[4] = 22./3;   d[4] = -16.;    e[4] = 32./3;
}

/*------------------------------------------------------------------*/

void BaseQ43d::point(const Vertex3d& s)const
{
  for(int i=0;i<NDOF;i++)
    {
      int ix = (i%NDOF2d) % NDOF1d;
      int iy = (i%NDOF2d) / NDOF1d;
      int iz = i/NDOF2d;

      N  [i]     = psi   (ix,s.x()) * psi   (iy,s.y()) * psi   (iz,s.z());
      DN [i].x() = psi_x (ix,s.x()) * psi   (iy,s.y()) * psi   (iz,s.z());
      DN [i].y() = psi   (ix,s.x()) * psi_x (iy,s.y()) * psi   (iz,s.z());
      DN [i].z() = psi   (ix,s.x()) * psi   (iy,s.y()) * psi_x (iz,s.z());
    }
}

/*------------------------------------------------------------------*/
}
#undef NDOF
#undef NDOF1d
#undef NDOF2d

