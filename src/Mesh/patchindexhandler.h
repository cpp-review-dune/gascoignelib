#ifndef  __PatchIndexHandler_h
#define  __PatchIndexHandler_h

#include  "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class PatchIndexHandler
{
protected:

  bool                               haspatch;
  nvector<IntVector>      indexofpatch;
  nvector<IntVector>      patch2cell;
  
  int dim;

public:

  int&                           GetDim()      { return dim;}
  bool&                          GetHasPatch() { return haspatch;}
  nvector<IntVector>& GetIndex()    { return indexofpatch;}

  IntVector&          GetPatch2Cell(int i)
    { assert(i<patch2cell.size()); return patch2cell[i];}
  const IntVector&    GetPatch2Cell(int i) const
    { assert(i<patch2cell.size()); return patch2cell[i];}
  
        nvector<IntVector>& GetAllPatch2Cell()       { return patch2cell; }
  const nvector<IntVector>& GetAllPatch2Cell() const { return patch2cell; }
  
  int npatches()    const { return indexofpatch.size();}
  bool HasPatch()   const { return haspatch;}

  const IntVector& IndicesOfPatch(int i) const { return indexofpatch[i];}
  IntVector CoarseIndices(int iq) const;

  int nodes_per_patch() const
    { 
      if (dim==2) return 9;
      return 27;
    }
};
}

#endif
