#ifndef __coarsehierarchicalmesh3d_h
#define __coarsehierarchicalmesh3d_h

#include  "hierarchicalmesh3d.h"

/*---------------------------------------------------*/

class CoarseHierarchicalMesh3d : public HierarchicalMesh3d
{
  protected :
    
  IntSet  CellRefList, CellCoarseList;
  IntVector   cn2o;

  void loop(IntVector& dst);

  public:
  
  CoarseHierarchicalMesh3d(const HierarchicalMesh3d&);
  void BasicInit(int depth);
  void GetRefinedList(IntVector&);
  void GetCoarsedList(IntVector&);
  void refine(const IntVector& cell_ref, const IntVector& cell_coarse);
};

/*---------------------------------------------------*/

#endif
