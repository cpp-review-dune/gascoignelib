#include  "gascoignemesh2d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
GascoigneMesh2d::GascoigneMesh2d()
{
}

/*-----------------------------------------*/

GascoigneMesh2d::~GascoigneMesh2d() 
{
}

/*---------------------------------------------------*/

IntVector GascoigneMesh2d::IndicesOfCell(int iq) const
{
  IntVector indices(4);

  int iq4 = iq*4;

  indices[0] = nc[iq4];
  indices[1] = nc[iq4+1];
  indices[2] = nc[iq4+3];
  indices[3] = nc[iq4+2];

  return indices;
}
}
