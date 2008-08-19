#ifndef  __hanglist_h
#define  __hanglist_h

#include  <string>
#include  "edgearray.h"
#include  "hang.h"

#include  <map>
#ifdef __OLDCOMPILER__
#include  <hash_map>
#define HANGMAP  hash_map<EdgeArray<N>,Hang,EdgeHash>
#else
#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HANGMAP   std::tr1::unordered_map<EdgeArray<N>,Hang,EdgeHash> 
#else
#include  <ext/hash_map>
#define HANGMAP  __gnu_cxx::hash_map<EdgeArray<N>,Hang,EdgeHash> 
#endif
#endif

/*------------------------------------------------------*/

namespace Gascoigne
{
class FixArrayHash
{
 public:
  template<int N>
    int operator()(const fixarray<N,int>& h) const { return h[0];}
};

/*------------------------------------------------------*/

//
/// This hash function has to be consistent with the operator "=="
/// for EdgeArrays, i.e. permutated fixarrays
//

class EdgeHash
{
 public:
  template<int N>
    int operator()(const EdgeArray<N>& h) const { return h.sum();}
};


/*------------------------------------------------------*/

template<int N>
/* class HangList : public std::map<EdgeArray<N>,Hang,EdgeArrayCompare<N> > */
class HangList : public HANGMAP
{
 protected:

 public:

  typedef typename HANGMAP::iterator        iterator;
  typedef typename HANGMAP::const_iterator  const_iterator;

  void update(const std::vector<int>&);
  void update(const std::vector<int>&, const std::vector<int>&);
  void make_consistent(HangList<N>&);
  void move(HangList<N>& src, iterator& p);
  HangList<N>& operator=(const HangList<N>& A);
  void BinWrite(std::ostream &s) const;
  void BinRead(std::istream &s);
};

/*------------------------------------------------------*/

template<int N>
std::ostream& operator<<(std::ostream &s, const HangList<N>& A);

template<int N>
std::istream& operator>>(std::istream &s, HangList<N>& A);
}

/*------------------------------------------------------*/

#undef HANGMAP

#endif
