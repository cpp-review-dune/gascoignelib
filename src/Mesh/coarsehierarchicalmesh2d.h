#ifndef __coarsehierarchicalmesh2d_h
#define __coarsehierarchicalmesh2d_h

#include  "hierarchicalmesh2d.h"

/*---------------------------------------------------*/

class CoarseHierarchicalMesh2d : public HierarchicalMesh2d
{
  protected :
    
  Gascoigne::IntSet  CellRefList, CellCoarseList;
  Gascoigne::IntVector   cn2o;

  void loop(Gascoigne::IntVector& dst);

  public:
  
  CoarseHierarchicalMesh2d(const HierarchicalMesh2d&);
  void BasicInit(int pdepth);
  void GetRefinedList(Gascoigne::IntVector&);
  void GetCoarsedList(Gascoigne::IntVector&);
  void refine(const Gascoigne::IntVector& cell_ref, const Gascoigne::IntVector& cell_coarse);
};

/*---------------------------------------------------*/

#endif
