#include "baseq23d.h"

#define NDOF   27
#define NDOF1d 3
#define NDOF2d 9

namespace Gascoigne
{

/*------------------------------------------------------------------*/

BaseQ23d::BaseQ23d()
{
  N  .resize(NDOF);  
  DN .resize(NDOF);  
  //DDN.resize(NDOF);

  //for(int i=0;i<NDOF;i++) DDN[i].zero();

  a[0] = 1.;  b[0] = -3.;  c[0] = 2.;
  a[1] = 0.;  b[1] =  4.;  c[1] = -4.;
  a[2] = 0.;  b[2] = -1.;  c[2] = 2.;
}

/*------------------------------------------------------------------*/

void BaseQ23d::point(const Vertex3d& s)const
{
  for(int i=0;i<NDOF;i++) 
    { 
      // nicht geprueft !!!!!!!!!!
//       int ix = i%NDOF1d;
//       int iy = i/NDOF1d;
//       int iz = (i%NDOF2d)/NDOF1d;

      int ix = (i%NDOF2d) % NDOF1d;
      int iy = (i%NDOF2d) / NDOF1d;
      int iz = i/NDOF2d;

      double px = psi(ix,s.x());
      double py = psi(iy,s.y());
      double pz = psi(iz,s.z());
      
      double dpx = psi_x(ix,s.x());
      double dpy = psi_x(iy,s.y());
      double dpz = psi_x(iz,s.z());
      
      N  [i]     =  px *  py *  pz;
      DN [i].x() = dpx *  py *  pz;
      DN [i].y() =  px * dpy *  pz;
      DN [i].z() =  px *  py * dpz;
      
//       double ddpx = psi_xx(ix,s.x());
//       double ddpy = psi_xx(iy,s.y());
//       double ddpz = psi_xx(iz,s.z());
      
//       DDN[i][0]  = ddpx *   py *   pz;
//       DDN[i][1]  =   px * ddpy *   pz;
//       DDN[i][2]  =   px *   py * ddpz;
      
//       DDN[i][3]  = dpx * dpy *  pz;
//       DDN[i][4]  = dpx *  py * dpz;
//       DDN[i][5]  =  px * dpy * dpz;
    }
}

/*------------------------------------------------------------------*/
}
#undef NDOF
#undef NDOF1d
#undef NDOF2d

