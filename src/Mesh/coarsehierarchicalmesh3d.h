#ifndef __coarsehierarchicalmesh3d_h
#define __coarsehierarchicalmesh3d_h

#include  "hierarchicalmesh3d.h"

/*---------------------------------------------------*/

class CoarseHierarchicalMesh3d : public HierarchicalMesh3d
{
  protected :
    
  Gascoigne::IntSet  CellRefList, CellCoarseList;
  Gascoigne::IntVector   cn2o;

  void loop(Gascoigne::IntVector& dst);

  public:
  
  CoarseHierarchicalMesh3d(const HierarchicalMesh3d&);
  void BasicInit(int depth);
  void GetRefinedList(Gascoigne::IntVector&);
  void GetCoarsedList(Gascoigne::IntVector&);
  void refine(const Gascoigne::IntVector& cell_ref, const Gascoigne::IntVector& cell_coarse);
};

/*---------------------------------------------------*/

#endif
