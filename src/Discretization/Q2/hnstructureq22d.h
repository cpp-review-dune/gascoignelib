#ifndef  __HNStructureQ22d_h
#define  __HNStructureQ22d_h

#include  "hnstructureq12d.h"

namespace Gascoigne
{

/*-----------------------------------------*/

class HNStructureQ22d : public  HNStructureQ12d
{
public:

  HNStructureQ22d();
    
  void Average   (GlobalVector& u) const;
  void Distribute(GlobalVector& u) const;
  void CondenseHanging(EntryMatrix& E, IntVector& indices) const;
  void CondenseHanging(IntVector& indices) const;
  void CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const { assert(0);}
};
}
#endif
