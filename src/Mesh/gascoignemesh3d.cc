#include  "gascoignemesh3d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
GascoigneMesh3d::GascoigneMesh3d()
{
}

/*---------------------------------------------------*/

IntVector GascoigneMesh3d::IndicesOfCell(int iq) const
{
  IntVector indices(8);
  
  int offset = 8*iq;

  indices[0] = nc[offset];
  indices[1] = nc[offset+1];
  indices[2] = nc[offset+3];
  indices[3] = nc[offset+2];
  indices[4] = nc[offset+4];
  indices[5] = nc[offset+5];
  indices[6] = nc[offset+7];
  indices[7] = nc[offset+6];

  return indices;
}
}
