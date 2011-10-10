#ifndef  __EdgeArray_h
#define  __EdgeArray_h

#include  "fixarray.h"

/*-----------------------------------------*/

namespace Gascoigne
{
template<int N>
class EdgeArray : public fixarray<N,int>
{
 public:

  EdgeArray<N>(const fixarray<N,int> &e);
  
  bool operator==(const fixarray<N,int> &e) const;

  int sum() const;
};

}

/*------------------------------------------------------*/

#endif
