#include "hanglist.h"
#include <cassert>

#ifdef __OLDCOMPILER__
#define HANGMAP  hash_map<EdgeArray<N>,Hang,EdgeHash>
#else
#define HANGMAP  __gnu_cxx::hash_map<EdgeArray<N>,Hang,EdgeHash> 
#endif

/*********************************************************************/

template<int N>
std::ostream& operator<<(std::ostream &s, const HangList<N>& A)
{
  s << A.size() << " hangs" << std::endl;
  typename HANGMAP::const_iterator p = A.begin();
  
  for (typename HANGMAP::const_iterator
	 p = A.begin(); p!=A.end(); p++)
    {
      s << p->first << "-> " << p->second;
    }
  return s;
}

/*********************************************************************/

template<int N>
std::istream& operator>>(std::istream &s, HangList<N>& A)
{
  int n;
  std::string symbol;
  
  s >> n >> symbol;

  assert(symbol=="hangs");

  fixarray<N,int>  ev;
  Hang         info;
  for (int i=0; i<n; i++)
    {
      s >> ev >> symbol >> info;
      A.insert(std::make_pair(ev,info));
    }
}

/*********************************************************************/

template std::ostream& operator<<(std::ostream &s, const HangList<2>& A);
template std::ostream& operator<<(std::ostream &s, const HangList<4>& A);
template std::istream& operator>>(std::istream &s, HangList<2>& A);
template std::istream& operator>>(std::istream &s, HangList<4>& A);

#undef HANGMAP
