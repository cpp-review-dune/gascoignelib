/*----------------------------   mesharray.h     ---------------------------*/
/*      $Id: mesharray.h,v 1.2 2008/09/05 16:04:25 richter Exp $                 */
#ifndef __mesharray_H
#define __mesharray_H
/*----------------------------   mesharray.h     ---------------------------*/


#include  <cstdlib> 
#include  <iostream> 
#include  <iterator> 
#include  <algorithm>
#include  <iostream>
#include  <cassert> 

/*-------------------------------------------------*/

namespace Tsuchimikado
{

  /**
   * array of fixed size with some usefull functions
   **/
  template<int N, class T>
    class mesharray 
    {
    public:
      
      typedef  T*        iterator;
      typedef  const T*  const_iterator;
  
    private:
      T __data[N];

    public:

      mesharray<N,T>()           { *this = T(); }
      mesharray<N,T>(const T& d) { *this = d; }
      
      const T*  begin() const { return &(__data[0]);}
      const T*  end  () const { return &(__data[0])+N;}
      T*        begin()       { return &(__data[0]);}
      T*        end  ()       { return &(__data[0])+N;}
      
      
      inline       T& operator[](int i)       { return __data[i]; }
      inline const T& operator[](int i) const { return __data[i]; }

      inline       T& operator()(int i)       { assert(i<N); return __data[i]; }
      inline const T& operator()(int i) const { assert(i<N); return __data[i]; }
      
      inline mesharray<N,T>& operator=(const T& d) 
	{
	  iterator  p(end());
	  while(p>begin()) *--p = d;
	  return *this;
	}

      inline mesharray<N,T>&  operator+=(const mesharray<N,T>& v)
	{ add(1.,v); return *this; }

      inline mesharray<N,T>&  operator*=(const T& d)
	{ iterator f = begin();
	  const_iterator l = end();
	  while (f!=l)
	    (*f++) *= d;
	  return *this;	}

      inline double  operator*(const mesharray<N,T>& d)
	{
	  double res=0;
	  iterator f1       = begin();
	  const_iterator f2 = d.begin();
	  const_iterator l1 = end();
	  while (f1!=l1)
	    res += (*f1++) * (*f2++);
	  return res;
	}

      inline void add(T d, const mesharray<N,T>& v)
      {
	iterator first  = mesharray<N,T>::begin();
	const_iterator last   = mesharray<N,T>::end();
	const_iterator first2 = v.mesharray<N,T>::begin();
	while(first != last)
	  {
	    (*first++) += d*(*first2++);
	  }
      }

      
/*       mesharray<N,T>& operator=(const fixarray<N,T>& v)  */
/* 	{ */
/* 	  iterator        p(begin()); */
/* 	  const_iterator  q(v.begin()); */
/* 	  while(p<end()) *p++ = *q++; */
/* 	  return *this; */
/* 	}  */


    };

  // some tools
  mesharray<9,int> getmesharray(int i1,int i2,int i3,int i4,int i5,int i6,int i7,int i8,int i9);
  mesharray<8,int> getmesharray(int i1,int i2,int i3,int i4,int i5,int i6,int i7,int i8);
  mesharray<6,int> getmesharray(int i1,int i2,int i3,int i4,int i5,int i6);
  mesharray<4,int> getmesharray(int i1,int i2,int i3,int i4);
  mesharray<3,int> getmesharray(int i1,int i2,int i3);
  mesharray<2,int> getmesharray(int i1,int i2);
  template<int N,class T> std::ostream& operator<<(std::ostream &s, const mesharray<N,T>& A);
  template<int N,class T> std::istream& operator>>(std::istream &s, mesharray<N,T>& A);
  
 
  
}


/*----------------------------   mesharray.h     ---------------------------*/
/* end of #ifndef __mesharray_H */
#endif
/*----------------------------   mesharray.h     ---------------------------*/
