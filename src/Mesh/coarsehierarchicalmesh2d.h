#ifndef __coarsehierarchicalmesh2d_h
#define __coarsehierarchicalmesh2d_h

#include  "hierarchicalmesh2d.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class CoarseHierarchicalMesh2d : public HierarchicalMesh2d
{
  protected :
    
  IntSet  CellRefList, CellCoarseList;
  IntVector   cn2o;

  void loop(IntVector& dst);

  public:
  
  CoarseHierarchicalMesh2d(const HierarchicalMesh2d&);
  void BasicInit(int pdepth);
  void GetRefinedList(IntVector&);
  void GetCoarsedList(IntVector&);
  void refine(const IntVector& cell_ref, const IntVector& cell_coarse);
};
}

/*---------------------------------------------------*/

#endif
