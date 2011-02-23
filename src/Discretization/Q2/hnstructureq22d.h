#ifndef  __HNStructureQ22d_h
#define  __HNStructureQ22d_h

#include  "hnstructureq12d.h"

namespace Gascoigne
{

/*-----------------------------------------*/

class HNStructureQ22d : public  HNStructureQ12d
{
protected:
  DoubleVector q1wei;

public:

  HNStructureQ22d();
    
  void Average   (GlobalVector& u) const;
  void Distribute(GlobalVector& u) const;
  void CondenseHanging(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingLowerHigher(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingHigherLower(EntryMatrix& E, IntVector& indices) const;
  void CondenseHangingMixed(EntryMatrix& E, IntVector& indices, int k) const;
  void CondenseHanging(IntVector& indices) const;
  void CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const { 
    std::cerr << "\"HNStructureQ22d::CondenseHangingPatch\" not written!" << std::endl;
    abort();
  }

  //void NewCondenseHanging(EntryMatrix& E, IntVector& indices1, IntVector& indices2) const;
};
}
#endif
