#ifndef  __compareclass_h
#define  __compareclass_h

#include "nvector.h"
#include <algorithm>

/*********************************************************/

template<class T>
class CompareLess
{
  protected:
  
  const T*   wp;

  public:

  CompareLess<T>() {}
  CompareLess<T>(const T& w)
    { 
      wp  = &w;
    }
  bool operator()(int i1,int i2) const
    {
      if( (*wp)[i1]  < (*wp)[i2] ) return 1;
      return 0;
    }
};

/*********************************************************/

template<class T>
class CompareObject
{
  protected:
  
  const T& P;

  public:

  CompareObject(const T& p) : P(p)
    {};
  
  bool operator()(int i, int j) const
    { return P[i]<P[j]; }
};

/*********************************************************/

template<class T>
class CompareObjectBigToSmall
{
  protected:
  
  const T& P;
  
  public:
  
  CompareObjectBigToSmall(const T& p) : P(p)
    {};
  
  bool operator()(int i, int j) const
    { return P[j]<P[i]; }
};

#endif

