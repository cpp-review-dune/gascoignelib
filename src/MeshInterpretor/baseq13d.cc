#include "baseq13d.h"

/*------------------------------------------------------------------*/

BaseQ13d::BaseQ13d() 
{
  BasicInit();
}

/*------------------------------------------------------------------*/

void BaseQ13d::point(const Vertex3d& s) const
{
  for(int i=0;i<8;i++) 
    { 
      int ix = i%2;
      int iy = (i%4)/2;
      int iz = i/4;

      double px = psi(ix,s.x());
      double py = psi(iy,s.y());
      double pz = psi(iz,s.z());

      double dx = psi_x(ix,s.x());
      double dy = psi_x(iy,s.y());
      double dz = psi_x(iz,s.z());

      N  [i]     = px * py * pz;
      DN [i].x() = dx * py * pz;
      DN [i].y() = px * dy * pz;
      DN [i].z() = px * py * dz;

//       dxy[i]     = dx * dy * pz;
//       dxz[i]     = dx * py * dz;
//       dyz[i]     = px * dy * dz;
//       dxyz[i]    = dx * dy * dz;
    }
}

/*------------------------------------------------------------------*/

