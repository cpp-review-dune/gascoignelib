#ifndef  __quadraticstructure_h
#define  __quadraticstructure_h

#include  <map>
#include  "fixarray.h"

/*--------------------------------------------------------------*/

template<int N>
class QuadraticHNStructure : public std::map<int,fixarray<N,int> >
{
};

#endif
