#ifndef  __vecalgo_h
#define  __vecalgo_h

#include  <vector>
#include  <set>
#include  "fadamath.h"

namespace Gascoigne
{
void transfer(int n, std::vector<int>& tr, const std::set<int>& del);
void transfer(int n, std::vector<int>& tr, std::vector  <int>& del);

template<class C>
void compress(std::vector<C>& dst, const std::vector<int>& src) {
  //int n = 0;
  int mmax = 0;

  for (int i=0; i<src.size(); i++)
    {
      int j = src[i];
      if (j>=0)
	{
	  dst[j] = dst[i];
	  mmax = max_int(mmax,j);
	  //n++;
	}
    }
  //dst.resize(n);
  dst.resize(mmax+1);
}
}

#endif
