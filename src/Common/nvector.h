#ifndef __nvector_h
#define __nvector_h

#include  <vector>
#include  <iterator>
#include  <iostream>
#include  <iterator>
#include  <climits>
#include  <numeric>
#include  <cassert>

#include  "fadamath.h"

/*----------------------------------------------*/

template<class T>
class nvector : public std::vector<T>
{
private:

public:

  typedef typename std::vector<T>::iterator       iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  ~nvector()  {}
  nvector()                    : std::vector<T>()    {}
  nvector(size_t n)               : std::vector<T>(n)   {}
  nvector(size_t n, const T& d)   : std::vector<T>(n,d) {}
  nvector(const std::vector<T>& v) : std::vector<T>(v)   {}

  friend std::ostream& operator<<(std::ostream &s, const nvector<T>& A)
  {
    std::ostream_iterator<T>  os(s, " ");
    copy(A.begin(),A.end(),os);
    return s;
  }
  friend std::istream& operator>>(std::istream &s, nvector<T>& A)
    {
      iterator p = A.begin();
      while(p!=A.end())
	s >> *p++;
      return s;
    }

  void write_data(std::ostream& s) const
  {
    s << size() << std::endl << *this;
  }

  void read_data(std::istream& s)
  {
    size_t n;
    s >> n;
    reservesize(n);
    s >> *this;
  }

  const T& secure_access(int i) const {
    assert(i<size());
    assert(i>=0);
    return operator[](i);
  }
  T& secure_access(int i) {
    assert(i<size());
    assert(i>=0);
    return operator[](i);
  }


  T operator*(const nvector& v) const;
  nvector<T>&   operator=(const T&);
  nvector<T>&   operator=(const std::vector<T>&);
  nvector<T>&   operator*=(const T& d);
  nvector<T>&   operator+=(const T& d) { add(d); return *this; }
  nvector<T>&   operator+=(const nvector& v) { add(1,v); return *this; }
  nvector<T>&   operator-=(const nvector& v) { add(-1,v); return *this; }

  void   zero ();
  void   equ  (const T&);
  void   equ  (const T&, const nvector&);
  void   equ  (const T&, const nvector&, const T&, const nvector&);
  void   equ  (const T&, const nvector&, const T&, const nvector&, 
	       const T&, const nvector&);
  void   sequ (const T&, const T&, const nvector&);
  void   add  (const T&);
  void   add  (const T&, const nvector&);
  void   add  (const T&, const nvector&, const T&, const nvector&);
  void   add  (const T&, const nvector&, const T&, const nvector&, 
	       const T&, const nvector&);
  void   sadd (const T&, const T&, const nvector&);
  double max  ()     const;
  double min  ()     const;
  T      sum()     const;
  double      norm()     const;
  double      norm_l1()     const;
  double      norm_l8()     const;

  void reservesize(size_t n)
    {
      std::vector<T>::reserve(n); std::vector<T>::resize(n);
    }
  void memory(size_t n)
    {
      std::vector<T>::reserve(n); std::vector<T>::resize(n);
    }
  void reservesize(size_t n, const T& s)
    {
      std::vector<T>::reserve(n); std::vector<T>::resize(n,s);
    }
  void reservesize(const nvector<T>& v) {reservesize(v.size());}

  void BinWrite(std::ostream& out) const;
  void BinRead (std::istream& in);
  int  find(const T& x) const;
};


/**************************************************/

template<class T>
inline double nvector<T>::norm() const
{
  const_iterator first  = begin();
  const_iterator last   = end();

  T n(0);
  while( first != last)
    {
      n += ((*first)) * ((*first));
      first++;
    }
  return sqrt(static_cast<double>(n));
}

/**************************************************/

template<class T>
inline T nvector<T>::sum() const
{
  const_iterator first  = begin();
  const_iterator last   = end();

  T n(0);
  while( first != last)
    {
      n += (*first++);
    }
  return n;
}

/**************************************************/

template<class T>
inline double nvector<T>::norm_l1() const
{
  const_iterator first  = begin();
  const_iterator last   = end();

  double n(0);
  while( first != last)
    {
      n += fabs((*first++));
    }
  return n;
}

/**************************************************/

template<class T>
inline double nvector<T>::norm_l8() const
{
  const_iterator first  = begin();
  const_iterator last   = end();

  double n(0);
  while( first != last)
    {
      n = GascoigneMath::max(n,fabs(*first));
      first++;
    }
  return n;
}

/**************************************************/

template<class T>
inline nvector<T>& nvector<T>::operator=(const T& d)
{
  iterator  first  = begin();
  const_iterator last   = end();

  while( first != last)
    {
      *first++ = d;
    }
  return *this;
}

/**************************************************/

template<class T>
inline nvector<T>& nvector<T>::operator=(const std::vector<T>& v)
{
  assert(size()==v.size());
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator vfirst = v.begin();

  while( first != last)
    {
      *first++ = *vfirst++;
    }
  return *this;
}

/**************************************************/

template<class T>
inline nvector<T>& nvector<T>::operator*=(const T& d)
{
  iterator  first  = begin();
  const_iterator last   = end();

  while(first != last)
    {
      (*first++) *= d;
   }
  return *this;
}

/**************************************************/

template<class T>
inline void nvector<T>::zero()
{
  iterator  first  = begin();
  const_iterator last   = end();

  while( first != last)
    {
      *first++ = 0;
    }
}

/**************************************************/

template<class T>
inline double nvector<T>::max() const
{
  double d = 0;//std::numeric_limits<double>::min();
/*   double d = std::numeric_limits<double>::min(); */
  const_iterator first  = begin();
  const_iterator last   = end();

  while( first != last)
    {
      d = GascoigneMath::max( d, fabs((*first)));
      first++;
    }
  return d;
}

/**************************************************/

template<class T>
inline double nvector<T>::min() const
{
  double d = 100000.;//std::numeric_limits<double>::max();
/*   double d = std::numeric_limits<double>::max(); */
  const_iterator first  = begin();
  const_iterator last   = end();

  while( first != last)
    {
      d = GascoigneMath::min( d, fabs((*first)));
      first++;
    }
  return d;
}

/**************************************************/

template<class T>
inline T nvector<T>::operator* (const nvector<T>& v) const
{
  const_iterator first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();

  T d(0);
  while(first != last)
    {
      d += (*first++) * (*first2++);
    }
  return d;

  //return inner_product(first,last,first2,0.);
}

/**************************************************/

template<class T>
inline void nvector<T>::equ (const T& d)
{
  iterator  first  = begin();
  const_iterator last   = end();

  while(first != last)
    {
      (*first++) = d;
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::equ (const T& d, const nvector<T>& v)
{
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++);
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::equ (const T& d, const nvector<T>& v, 
			     const T& e, const nvector<T>& w)
{
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();
  const_iterator first3 = w.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++);
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::equ (const T& d, const nvector<T>& v, 
			     const T& e, const nvector<T>& w,
			     const T& f, const nvector<T>& x)
{
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();
  const_iterator first3 = w.begin();
  const_iterator first4 = x.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++) + f*(*first4++);
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::sequ (const T& s, const T& d, const nvector<T>& v)
{
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();

  while(first != last)
    {
      (*first) = s*(*first) + d*(*first2++);
      first++;
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::add (const T& d)
{
  iterator  first  = begin();
  const_iterator last   = end();

  while(first != last)
    {
      (*first++) += d;
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::add (const T& d, const nvector<T>& v)
{
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++);
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::add (const T& d, const nvector<T>& v,
			     const T& e, const nvector<T>& w)
{
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();
  const_iterator first3 = w.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++) + e*(*first3++);
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::add (const T& d, const nvector<T>& v,
			     const T& e, const nvector<T>& w,
			     const T& f, const nvector<T>& x)
{
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();
  const_iterator first3 = w.begin();
  const_iterator first4 = x.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++) + e*(*first3++) + f*(*first4++);
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::sadd (const T& a, const T& d, const nvector<T>& v)
{
  iterator  first  = begin();
  const_iterator last   = end();
  const_iterator first2 = v.begin();

  while(first != last)
    {
      (*first) = a*(*first) + d*(*first2++);
      first++;
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::BinWrite(std::ostream& out) const
{
//#ifdef __OLDCOMPILER__
//  out << size() << std::endl << "[";
//  int laenge = reinterpret_cast<const char*>(end()) - reinterpret_cast<const char*>(begin());
//  out.write (reinterpret_cast<const char*>(begin()),laenge);
//  out << "]";  
//#else
//  std::cerr << "nvector<T>::BinWrite\n\tkaputt wegen gcc 3.1\n"; exit(1);
//#endif
  out << size() << std::endl << "[";
  for(int i=0; i<size(); i++)
    {
      out.write (reinterpret_cast<const char*>(&(operator[](i))),sizeof(operator[](i)));
    }
  out << "]"; 
}

/**********************************************************/

template<class T>
inline void nvector<T>::BinRead(std::istream& in)
{
//#ifdef __OLDCOMPILER__
//  char c;
//  int  n;
//  in >> n >> c;
//  resize(n);
//  int laenge = reinterpret_cast<const char*>(end()) - reinterpret_cast<const char*>(begin());
//  in.read (reinterpret_cast<void*>(begin()),laenge);
//  in >> c;  
//#else
//  std::cerr << "nvector<T>::BinRead\n\tkaputt wegen gcc 3.1\n"; exit(1);
//#endif
  char c;
  int  n;
  in >> n >> c;
  resize(n);
  for(int i=0; i<size(); i++)
    {
      in.read(reinterpret_cast<char*>(&(operator[](i))),sizeof(operator[](i)));
    }
  in >> c;
}

/**********************************************************/

template<class T>
inline int nvector<T>::find(const T& x) const
{
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i]==x) return i;
    }
  return -1;
}

/**********************************************************/

#endif

