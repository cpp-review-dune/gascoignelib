#include  "patchindexhandler.h"

/*-----------------------------------------*/

nvector<int> PatchIndexHandler::CoarseIndices(int iq) const
{
  nvector<int> indices;

  assert(iq>=0);
  assert(iq<indexofpatch.size());

  if (dim==2)
    {
      indices.resize(4);

      indices[0] = indexofpatch[iq][0];
      indices[1] = indexofpatch[iq][2];
      indices[2] = indexofpatch[iq][6];
      indices[3] = indexofpatch[iq][8];
    }
  else
    {
      indices.resize(8);

      indices[0] = indexofpatch[iq][0];
      indices[1] = indexofpatch[iq][2];
      indices[2] = indexofpatch[iq][6];
      indices[3] = indexofpatch[iq][8];
      indices[4] = indexofpatch[iq][18];
      indices[5] = indexofpatch[iq][20];
      indices[6] = indexofpatch[iq][24];
      indices[7] = indexofpatch[iq][26];
    }
  return indices;
}
