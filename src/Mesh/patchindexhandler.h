#ifndef  __PatchIndexHandler_h
#define  __PatchIndexHandler_h

#include  "gascoigne.h"

/*-----------------------------------------*/

class PatchIndexHandler
{
protected:

  bool                               haspatch;
  nvector<Gascoigne::IntVector>      indexofpatch;
  nvector<Gascoigne::IntVector>      patch2cell;
  
  int dim;

public:

  int&                           GetDim()      { return dim;}
  bool&                          GetHasPatch() { return haspatch;}
  nvector<Gascoigne::IntVector>& GetIndex()    { return indexofpatch;}

  Gascoigne::IntVector&          GetPatch2Cell(int i)
    { assert(i<patch2cell.size()); return patch2cell[i];}
  const Gascoigne::IntVector&    GetPatch2Cell(int i) const
    { assert(i<patch2cell.size()); return patch2cell[i];}
  
        nvector<Gascoigne::IntVector>& GetAllPatch2Cell()       { return patch2cell; }
  const nvector<Gascoigne::IntVector>& GetAllPatch2Cell() const { return patch2cell; }
  
  int npatches()    const { return indexofpatch.size();}
  bool HasPatch()   const { return haspatch;}

  const nvector<int>& IndicesOfPatch(int i) const { return indexofpatch[i];}
  nvector<int> CoarseIndices(int iq) const;

  int nodes_per_patch() const
    { 
      if (dim==2) return 9;
      return 27;
    }
};

#endif
