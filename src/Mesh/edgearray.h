#ifndef  __EdgeArray_h
#define  __EdgeArray_h

#include  <map>
#include  "fixarray.h"

/*-----------------------------------------*/

namespace Gascoigne
{
template<int N>
class EdgeArray : public fixarray<N,int>
{
 public:

  EdgeArray<N>();
  EdgeArray<N>(const int& d);
  EdgeArray<N>(const fixarray<N,int>& e);
  EdgeArray<N>(const EdgeArray<N>& e);
  
  bool operator==(const EdgeArray<N>&)    const;
  bool operator==(const fixarray<N,int>&) const;

  int sum() const;
};

/*-----------------------------------------*/

template<int N>
class EdgeArrayCompare
{
 public:

  bool operator()(const EdgeArray<N>& EA1, const EdgeArray<N>& EA2) const
  {
/*     EdgeArray<N> S1(EA1), S2(EA2); */
/*     std::sort(S1.begin(),S1.end()); */
/*     std::sort(S2.begin(),S2.end()); */
/*     return S1<S2; */

    return !(EA1==EA2);
  }
};
}

/*------------------------------------------------------*/

#endif
