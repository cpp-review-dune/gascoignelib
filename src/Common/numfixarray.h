#ifndef __numfixarray_h
#define __numfixarray_h

#include  "fixarray.h"
#include  "fadamath.h"
#include  "nvector.h"

/*-------------------------------------------------*/

template<int N, class T>
class numfixarray  : public fixarray<N,T>
{
    typedef typename fixarray<N,T>::iterator         nvp;
    typedef typename fixarray<N,T>::const_iterator   cnvp;

  public:
    ~numfixarray()    {}
    numfixarray()                    : fixarray<N,T>(0.)    {}
/*   numfixarray(size_t n)               : fixarray<N,T>(n)   {} */
    numfixarray(const T& d)   : fixarray<N,T>(d) {}
    numfixarray(const numfixarray<N,T>& v) : fixarray<N,T>(v)   {}

    double operator*(const numfixarray& v) const;
    double operator*(const nvector<T>& v) const;
    numfixarray<N,T>&   operator=(const T&);
    numfixarray<N,T>&   operator=(const numfixarray&);
    numfixarray<N,T>&   operator*=(const numfixarray&);
    numfixarray<N,T>&   operator*=(double d);
    numfixarray<N,T>&   operator/=(double d);
    numfixarray<N,T>&   operator+=(double d) { add(d); return *this; }
    numfixarray<N,T>&   operator+=(const numfixarray& v) { add(1.,v); return *this; }
    numfixarray<N,T>&   operator-=(const numfixarray& v) { add(-1.,v); return *this; }
    numfixarray<N,T>    operator-(numfixarray v1) const
      { v1.add(-1.,*this);
	v1*=(-1);
	return v1;
      }
    
    int     n    () const { return size(); }
    void    zero ();
    void    equ  (const T&);
    void    equ  (double, const numfixarray&);
    void    equ  (double, const numfixarray&, double, const numfixarray&);
    void    equ  (double, const numfixarray&, double, const numfixarray&, 
		  double, const numfixarray&);
    void    equ  (double, const numfixarray&, double, const numfixarray&, 
		  double, const numfixarray&, double, const numfixarray&);
    void    equ  (double, const numfixarray&, double, const numfixarray&, 
		  double, const numfixarray&, double, const numfixarray&,
		  double, const numfixarray&, double, const numfixarray&);
    void    sequ (double, double, const numfixarray&);
    void    add  (double);
    void    add  (T, const numfixarray&);
    void    add  (double, const numfixarray&, double, const numfixarray&);
    void    add  (double, const numfixarray&, double, const numfixarray&, 
		  double, const numfixarray&);
    void    sadd (double, double, const numfixarray&);
    double  max  ()     const;
    double  norm ()     const;
    double  norm_l8 ()     const;
    double  norm_l1()     const;
    double  norm_l2()     const {return norm();}

    void normalise()
      {
	T d=0.;
	for(int i=0;i<n();i++)
	  {
	    d += (*this)[i]*(*this)[i];
	  }
	d = 1./sqrt(d);
	for(int i=0;i<n();i++)
	  {
	    (*this)[i] *= d;
	  }
      }

};


/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::norm_l8() const
{
  cnvp first  = begin();
  cnvp last   = end();

  double d=0.;
  while( first != last)
    {
      d = GascoigneMath::max(d,fabs(*first));
      first++;
    }
  return d;
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::norm() const
{
  cnvp first  = begin();
  cnvp last   = end();

  double n=0.;
  while( first != last)
    {
      n += ((*first)) * ((*first));
      first++;
    }
  return sqrt(n);
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::norm_l1() const
{
  cnvp first  = begin();
  cnvp last   = end();

  double n=0.;
  while( first != last)
    {
      n += fabs((*first++));
    }
  return n;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator=(const T& d)
{
  nvp  first  = begin();
  cnvp last   = end();

  while( first != last)
    {
      *first++ = d;
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator=(const numfixarray<N,T>& v)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp vfirst = v.begin();

  while( first != last)
    {
      *first++ = *vfirst++;
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator*=(const numfixarray<N,T>& d)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp fd     = d.begin();

  while(first != last)
    {
      (*first++) *= (*fd++);
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator*=(double d)
{
  nvp  first  = begin();
  cnvp last   = end();

  while(first != last)
    {
      (*first++) *= d;
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator/=(double d)
{
  nvp  first  = begin();
  cnvp last   = end();

  while(first != last)
    {
      (*first++) /= d;
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::zero()
{
  nvp  first  = begin();
  cnvp last   = end();

  while( first != last)
    {
      *first++ = 0.;
    }
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::max() const
{
  double d = -1.;
  cnvp first  = begin();
  cnvp last   = end();

  while( first != last)
    {
      d = MAX( d, fabs((*first)));
      first++;
    }
  return d;
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::operator* (const numfixarray<N,T>& v) const
{
  cnvp first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();

  double d = 0.;
  while(first != last)
    {
      d += (*first++) * (*first2++);
    }
  return d;
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::operator* (const nvector<T>& v) const
{
  cnvp first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();

  double d = 0.;
  while(first != last)
    {
      d += (*first++) * (*first2++);
    }
  return d;
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (const T& d)
{
  nvp  first  = begin();
  cnvp last   = end();

  while(first != last)
    {
      (*first++) = d;
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v, 
				   double e, const numfixarray<N,T>& w)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();
  cnvp first3 = w.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v, 
				   double e, const numfixarray<N,T>& w,
				   double f, const numfixarray<N,T>& x)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();
  cnvp first3 = w.begin();
  cnvp first4 = x.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++) + f*(*first4++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v, 
				   double e, const numfixarray<N,T>& w,
				   double f, const numfixarray<N,T>& x,
				   double g, const numfixarray<N,T>& y	   )
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();
  cnvp first3 = w.begin();
  cnvp first4 = x.begin();
  cnvp first5 = y.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++) + f*(*first4++) + g*(*first5++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v, 
				   double e, const numfixarray<N,T>& w,
				   double f, const numfixarray<N,T>& x,
				   double g, const numfixarray<N,T>& y,
				   double h, const numfixarray<N,T>& z,
				   double i, const numfixarray<N,T>& zz)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();
  cnvp first3 = w.begin();
  cnvp first4 = x.begin();
  cnvp first5 = y.begin();
  cnvp first6 = z.begin();
  cnvp first7 = zz.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++) + f*(*first4++) + g*(*first5++) 
		   + h*(*first5++) + i*(*first6++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::sequ (double s, double d, const numfixarray<N,T>& v)
{
  equ(s,*this);
  add(d,v);
  return;

  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();

  while(first != last)
    {
      (*first) = s*(*first) + d*(*first2++);
      first++;
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::add (double d)
{
  nvp  first  = begin();
  cnvp last   = end();

  while(first != last)
    {
      (*first++) += d;
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::add (T d, const numfixarray<N,T>& v)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::add (double d, const numfixarray<N,T>& v,
				   double e, const numfixarray<N,T>& w)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();
  cnvp first3 = w.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++) + e*(*first3++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::add (double d, const numfixarray<N,T>& v,
				   double e, const numfixarray<N,T>& w,
				   double f, const numfixarray<N,T>& x)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();
  cnvp first3 = w.begin();
  cnvp first4 = x.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++) + e*(*first3++) + f*(*first4++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::sadd (double a, double d, const numfixarray<N,T>& v)
{
  nvp  first  = begin();
  cnvp last   = end();
  cnvp first2 = v.begin();

  while(first != last)
    {
      (*first) = a*(*first) + d*(*first2++);
      first++;
    }
}

#endif

