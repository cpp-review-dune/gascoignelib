#ifndef  __HNStructureQ23d_h
#define  __HNStructureQ23d_h

#include  "hnstructureq13d.h"

namespace Gascoigne
{

/*-----------------------------------------*/

class HNStructureQ23d : public  HNStructureQ13d
{
  fixarray<9,double>  fwei;

  fixarray<12,fixarray<3,int> >  lnoe;
  fixarray< 6,fixarray<5,int> >  lnop;

  void CondenseHanging2er(EntryMatrix& E, nvector<int>& indices) const;
  void CondenseHanging4er(EntryMatrix& E, nvector<int>& indices) const;

public:

  HNStructureQ23d();

  void Average(GlobalVector& u) const;
  void Distribute(GlobalVector& u) const;
  void CondenseHanging(EntryMatrix& E, IntVector& indices) const;
  void CondenseHanging(IntVector& indices) const;
  void CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const { assert(0);}
};

}
#endif
