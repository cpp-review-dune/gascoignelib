#ifndef  __HNStructureQ1_h
#define  __HNStructureQ1_h

#include  "hnstructureinterface.h"
#include  <map>

/*-----------------------------------------*/

namespace Gascoigne
{
class HNStructureQ1 : public virtual HNStructureInterface
{
protected:

  typedef   fixarray<3,int>        EdgeVector;
  typedef   std::map<int,EdgeVector>::iterator        iterator;
  typedef   std::map<int,EdgeVector>::const_iterator  const_iterator;

public:

  HNStructureQ1() {}
  ~HNStructureQ1() {}
};
}

/*-----------------------------------------*/


#endif
