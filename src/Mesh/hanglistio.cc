#include "hanglist.h"
#include <cassert>

using namespace std;

#ifdef __OLDCOMPILER__
#define HANGMAP  hash_map<EdgeArray<N>,Hang,EdgeHash>
#else
#ifdef __NEWER_THAN_GCC_4_2__
#define HANGMAP  std::tr1::unordered_map<EdgeArray<N>,Hang,EdgeHash> 
#else
#define HANGMAP  __gnu_cxx::hash_map<EdgeArray<N>,Hang,EdgeHash> 
#endif
#endif

/*********************************************************************/

namespace Gascoigne
{
template<int N>
ostream& operator<<(ostream &s, const HangList<N>& A)
{
  s << A.size() << " hangs" << endl;
  typename HANGMAP::const_iterator p = A.begin();
  
  for (p = A.begin(); p!=A.end(); p++)
    {
      s << p->first << "-> " << p->second;
    }
  return s;
}

/*********************************************************************/

template<int N>
istream& operator>>(istream &s, HangList<N>& A)
{
  int n;
  string symbol;
  
  s >> n >> symbol;

  assert(symbol=="hangs");

  fixarray<N,int>  ev;
  Hang         info;
  for (int i=0; i<n; i++)
    {
      s >> ev >> symbol >> info;
      A.insert(make_pair(ev,info));
    }

  return s; 
}

/*********************************************************************/

template ostream& operator<<(ostream &s, const HangList<2>& A);
template ostream& operator<<(ostream &s, const HangList<4>& A);
template istream& operator>>(istream &s, HangList<2>& A);
template istream& operator>>(istream &s, HangList<4>& A);
}

#undef HANGMAP
