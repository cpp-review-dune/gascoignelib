#ifndef __find_in_linehang_h
#define __find_in_linehang_h

/*---------------------------------------------------*/

namespace Gascoigne
{
template <int N>
std::pair<typename HangList<N>::iterator,bool> 
find_in_linehang(HangList<N>& LineHang, const fixarray<N,int>& lineglob)
{
  // sort(lineglob.begin(),lineglob.end());
  typename HangList<N>::iterator p = LineHang.find(lineglob);
  bool b=0;
  if(p!=LineHang.end()) b = 1;
  return std::make_pair(p,b);
}
}

#endif
