#ifndef  __levelsorter_h
#define  __levelsorter_h

/*---------------------------------------------------*/

class LevelSorter2d
{
protected:

  const HierarchicalMesh2d& HM;

public:

  LevelSorter2d(const HierarchicalMesh2d& HMHM): HM(HMHM) {}
  bool operator()(int i, int j) const
    {
      return ( HM.quad(i).level() > HM.quad(j).level() );
    }
};

/*---------------------------------------------------*/

class HangEdgeSort3
{
protected:

  const LevelMesh2d& LR;

public:

  HangEdgeSort3(const LevelMesh2d& L) : LR(L) {}
  bool operator() (int i, int j) const
    {
      return (!LR.EdgeIsHangingGlobalIndex(i) && LR.EdgeIsHangingGlobalIndex(j));
    }
};

#endif
