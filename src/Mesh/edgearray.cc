#include  "edgearray.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{
template<int N>
EdgeArray<N>::EdgeArray<N>() : fixarray<N,int>() {}

template<int N>
EdgeArray<N>::EdgeArray<N>(const int& d) : fixarray<N,int>(d) {}

template<int N>
EdgeArray<N>::EdgeArray<N>(const fixarray<N,int>& e) : fixarray<N,int>(e) {};

template<int N>
EdgeArray<N>::EdgeArray<N>(const EdgeArray<N>& e) : fixarray<N,int>(e) {};

/*--------------------------------------------------------------*/

template class EdgeArray<2>;
template class EdgeArray<4>;

/*--------------------------------------------------------------*/

bool EdgeArray<2>::operator==(const EdgeArray<2>& A) const
{
  if ( ((*this)[0]!=A[0]) && ((*this)[0]!=A[1]) ) return 0;
  if ( ((*this)[1]!=A[0]) && ((*this)[1]!=A[1]) ) return 0;
  return 1;
}

/*--------------------------------------------------------------*/

bool EdgeArray<2>::operator==(const fixarray<2,int>& A) const
{
  if ( ((*this)[0]!=A[0]) && ((*this)[0]!=A[1]) ) return 0;
  if ( ((*this)[1]!=A[0]) && ((*this)[1]!=A[1]) ) return 0;
  return 1;
}

/*--------------------------------------------------------------*/

int EdgeArray<2>::sum() const 
{ 
  return (*this)[0]+(*this)[1];
}

/*--------------------------------------------------------------*/

bool EdgeArray<4>::operator==(const EdgeArray<4>& A) const
{
  if ( ((*this)[0]!=A[0]) && ((*this)[0]!=A[1]) && 
       ((*this)[0]!=A[2]) && ((*this)[0]!=A[3]) ) return 0;
  if ( ((*this)[1]!=A[0]) && ((*this)[1]!=A[1]) && 
       ((*this)[1]!=A[2]) && ((*this)[1]!=A[3]) ) return 0;
  if ( ((*this)[2]!=A[0]) && ((*this)[2]!=A[1]) && 
       ((*this)[2]!=A[2]) && ((*this)[2]!=A[3]) ) return 0;
  return 1;
}

/*--------------------------------------------------------------*/

bool EdgeArray<4>::operator==(const fixarray<4,int>& A) const
{
  if ( ((*this)[0]!=A[0]) && ((*this)[0]!=A[1]) && 
       ((*this)[0]!=A[2]) && ((*this)[0]!=A[3]) ) return 0;
  if ( ((*this)[1]!=A[0]) && ((*this)[1]!=A[1]) && 
       ((*this)[1]!=A[2]) && ((*this)[1]!=A[3]) ) return 0;
  if ( ((*this)[2]!=A[0]) && ((*this)[2]!=A[1]) && 
       ((*this)[2]!=A[2]) && ((*this)[2]!=A[3]) ) return 0;
  return 1;
}

/*--------------------------------------------------------------*/

int EdgeArray<4>::sum() const 
{ 
  return (*this)[0]+(*this)[1]+(*this)[2]+(*this)[3];
}
}

/*--------------------------------------------------------------*/

