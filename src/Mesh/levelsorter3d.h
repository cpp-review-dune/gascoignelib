#ifndef  __levelsorter3d_h
#define  __levelsorter3d_h

#include  "hierarchicalmesh3d.h"

/*---------------------------------------------------*/

class LevelSorter3d
{
protected:

  const HierarchicalMesh3d& HM;

public:

  LevelSorter3d(const HierarchicalMesh3d& HMHM): HM(HMHM) {}
  bool operator()(int i, int j) const
    {
      return ( HM.hex(i).level() > HM.hex(j).level() );
    }
};

/*---------------------------------------------------*/

#endif
