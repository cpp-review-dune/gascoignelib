#include  "baseq12dwith.h"

using namespace std;

#define NDOF   4
#define NDOF1d 2

/*------------------------------------------------------------------*/

BaseQ12dWith::BaseQ12dWith()  : BaseQ12d()
{
  dxy.reservesize(NDOF);
}

/*------------------------------------------------------------------*/

void BaseQ12dWith::point(const Vertex2d& s)
{
  BaseQ12d::point(s);
  for(int i=0;i<NDOF;i++) 
    { 
      int ix = i%NDOF1d;
      int iy = i/NDOF1d;

      dxy[i]     = psi_x(ix,s.x()) * psi_x(iy,s.y());
    }
}

/*------------------------------------------------------------------*/

#undef NDOF
#undef NDOF1d
