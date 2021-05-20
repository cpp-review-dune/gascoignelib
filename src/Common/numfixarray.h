/**
 *
 * Copyright (C) 2004, 2005 by the Gascoigne 3D authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/

#ifndef __numfixarray_h
#define __numfixarray_h

//#include  "fixarray.h"
#include "nvector.h"
#include <array>

/*-------------------------------------------------*/

namespace Gascoigne {
template<size_t N, typename T>
class numfixarray : public std::array<T, N>
{
  typedef typename std::array<T, N>::iterator nvp;
  typedef typename std::array<T, N>::const_iterator cnvp;

public:
  ~numfixarray() {}
  numfixarray() { std::array<T, N>::fill(0.0); }
  /*   numfixarray(size_t n)               : std::array<T,N>(n)   {} */
  numfixarray(const T& d) { std::array<T, N>::fill(d); }
  numfixarray(const numfixarray<N, T>& v)
    : std::array<T, N>(v)
  {}

  // Output
  template<size_t NN, typename TT>
  friend std::ostream& operator<<(std::ostream& os,
                                  const numfixarray<NN, TT>& X);

  double operator*(const numfixarray& v) const;
  double operator*(const nvector<T>& v) const;
  numfixarray<N, T>& operator=(const T&);
  numfixarray<N, T>& operator=(const numfixarray&);
  numfixarray<N, T>& operator*=(const numfixarray&);
  numfixarray<N, T>& operator*=(double d);
  numfixarray<N, T>& operator/=(double d);
  numfixarray<N, T>& operator+=(double d)
  {
    add(d);
    return *this;
  }
  numfixarray<N, T>& operator+=(const numfixarray& v)
  {
    add(1., v);
    return *this;
  }
  numfixarray<N, T>& operator-=(const numfixarray& v)
  {
    add(-1., v);
    return *this;
  }
  numfixarray<N, T> operator-(numfixarray v1) const
  {
    v1.add(-1., *this);
    v1 *= (-1);
    return v1;
  }

  int n() const { return std::array<T, N>::size(); }
  void zero();
  void equ(const T&);
  void equ(double, const numfixarray&);
  void equ(double, const numfixarray&, double, const numfixarray&);
  void equ(double,
           const numfixarray&,
           double,
           const numfixarray&,
           double,
           const numfixarray&);
  void equ(double,
           const numfixarray&,
           double,
           const numfixarray&,
           double,
           const numfixarray&,
           double,
           const numfixarray&);
  void equ(double,
           const numfixarray&,
           double,
           const numfixarray&,
           double,
           const numfixarray&,
           double,
           const numfixarray&,
           double,
           const numfixarray&,
           double,
           const numfixarray&);
  void equ(const std::array<double, 8>& s,
           const std::array<numfixarray<N, T>, 8>& x);
  void sequ(double, double, const numfixarray&);
  void add(double);
  void add(T, const numfixarray&);
  void add(double, const numfixarray&, double, const numfixarray&);
  void add(double,
           const numfixarray&,
           double,
           const numfixarray&,
           double,
           const numfixarray&);
  void sadd(double, double, const numfixarray&);
  double max() const;
  double norm() const;
  double norm_l8() const;
  double norm_l1() const;
  double norm_l2() const { return norm(); }

  void normalise()
  {
    T d = 0.;
    for (int i = 0; i < n(); i++) {
      d += (*this)[i] * (*this)[i];
    }
    d = 1. / sqrt(d);
    for (int i = 0; i < n(); i++) {
      (*this)[i] *= d;
    }
  }
};

template<size_t N, typename T>
std::ostream&
operator<<(std::ostream& os, const numfixarray<N, T>& X)
{
  for (auto it : X)
    os << it << "\t";
  return os;
}

/**************************************************/

template<size_t N, typename T>
inline double
numfixarray<N, T>::norm_l8() const
{
  cnvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  double d = 0.;
  while (first != last) {
    d = std::max(d, fabs(*first));
    first++;
  }
  return d;
}

/**************************************************/

template<size_t N, typename T>
inline double
numfixarray<N, T>::norm() const
{
  cnvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  double n = 0.;
  while (first != last) {
    n += ((*first)) * ((*first));
    first++;
  }
  return sqrt(n);
}

/**************************************************/

template<size_t N, typename T>
inline double
numfixarray<N, T>::norm_l1() const
{
  cnvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  double n = 0.;
  while (first != last) {
    n += fabs((*first++));
  }
  return n;
}

/**************************************************/

template<size_t N, typename T>
inline numfixarray<N, T>&
numfixarray<N, T>::operator=(const T& d)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  while (first != last) {
    *first++ = d;
  }
  return *this;
}

/**************************************************/

template<size_t N, typename T>
inline numfixarray<N, T>&
numfixarray<N, T>::operator=(const numfixarray<N, T>& v)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp vfirst = v.std::array<T, N>::begin();

  while (first != last) {
    *first++ = *vfirst++;
  }
  return *this;
}

/**************************************************/

template<size_t N, typename T>
inline numfixarray<N, T>&
numfixarray<N, T>::operator*=(const numfixarray<N, T>& d)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp fd = d.std::array<T, N>::begin();

  while (first != last) {
    (*first++) *= (*fd++);
  }
  return *this;
}

/**************************************************/

template<size_t N, typename T>
inline numfixarray<N, T>&
numfixarray<N, T>::operator*=(double d)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  while (first != last) {
    (*first++) *= d;
  }
  return *this;
}

/**************************************************/

template<size_t N, typename T>
inline numfixarray<N, T>&
numfixarray<N, T>::operator/=(double d)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  while (first != last) {
    (*first++) /= d;
  }
  return *this;
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::zero()
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  while (first != last) {
    *first++ = 0.;
  }
}

/**************************************************/

template<size_t N, typename T>
inline double
numfixarray<N, T>::max() const
{
  double d = -1.;
  cnvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  while (first != last) {
    d = std::max(d, fabs((*first)));
    first++;
  }
  return d;
}

/**************************************************/

template<size_t N, typename T>
inline double
numfixarray<N, T>::operator*(const numfixarray<N, T>& v) const
{
  cnvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();

  double d = 0.;
  while (first != last) {
    d += (*first++) * (*first2++);
  }
  return d;
}

/**************************************************/

template<size_t N, typename T>
inline double
numfixarray<N, T>::operator*(const nvector<T>& v) const
{
  cnvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.begin();

  double d = 0.;
  while (first != last) {
    d += (*first++) * (*first2++);
  }
  return d;
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::equ(const T& d)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  while (first != last) {
    (*first++) = d;
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::equ(double d, const numfixarray<N, T>& v)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();

  while (first != last) {
    (*first++) = d * (*first2++);
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::equ(double d,
                       const numfixarray<N, T>& v,
                       double e,
                       const numfixarray<N, T>& w)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();
  cnvp first3 = w.std::array<T, N>::begin();

  while (first != last) {
    (*first++) = d * (*first2++) + e * (*first3++);
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::equ(double d,
                       const numfixarray<N, T>& v,
                       double e,
                       const numfixarray<N, T>& w,
                       double f,
                       const numfixarray<N, T>& x)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();
  cnvp first3 = w.std::array<T, N>::begin();
  cnvp first4 = x.std::array<T, N>::begin();

  while (first != last) {
    (*first++) = d * (*first2++) + e * (*first3++) + f * (*first4++);
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::equ(double d,
                       const numfixarray<N, T>& v,
                       double e,
                       const numfixarray<N, T>& w,
                       double f,
                       const numfixarray<N, T>& x,
                       double g,
                       const numfixarray<N, T>& y)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();
  cnvp first3 = w.std::array<T, N>::begin();
  cnvp first4 = x.std::array<T, N>::begin();
  cnvp first5 = y.std::array<T, N>::begin();

  while (first != last) {
    (*first++) =
      d * (*first2++) + e * (*first3++) + f * (*first4++) + g * (*first5++);
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::equ(double d,
                       const numfixarray<N, T>& v,
                       double e,
                       const numfixarray<N, T>& w,
                       double f,
                       const numfixarray<N, T>& x,
                       double g,
                       const numfixarray<N, T>& y,
                       double h,
                       const numfixarray<N, T>& z,
                       double i,
                       const numfixarray<N, T>& zz)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();
  cnvp first3 = w.std::array<T, N>::begin();
  cnvp first4 = x.std::array<T, N>::begin();
  cnvp first5 = y.std::array<T, N>::begin();
  cnvp first6 = z.std::array<T, N>::begin();
  cnvp first7 = zz.std::array<T, N>::begin();

  while (first != last) {
    (*first++) = d * (*first2++) + e * (*first3++) + f * (*first4++) +
                 g * (*first5++) + h * (*first6++) + i * (*first7++);
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::equ(const std::array<double, 8>& s,
                       const std::array<numfixarray<N, T>, 8>& x)
{
  zero();
  for (int i = 0; i < 8; ++i)
    add(s[i], x[i]);
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::sequ(double s, double d, const numfixarray<N, T>& v)
{
  equ(s, *this);
  add(d, v);
  return;

  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();

  while (first != last) {
    (*first) = s * (*first) + d * (*first2++);
    first++;
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::add(double d)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();

  while (first != last) {
    (*first++) += d;
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::add(T d, const numfixarray<N, T>& v)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();

  while (first != last) {
    (*first++) += d * (*first2++);
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::add(double d,
                       const numfixarray<N, T>& v,
                       double e,
                       const numfixarray<N, T>& w)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();
  cnvp first3 = w.std::array<T, N>::begin();

  while (first != last) {
    (*first++) += d * (*first2++) + e * (*first3++);
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::add(double d,
                       const numfixarray<N, T>& v,
                       double e,
                       const numfixarray<N, T>& w,
                       double f,
                       const numfixarray<N, T>& x)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();
  cnvp first3 = w.std::array<T, N>::begin();
  cnvp first4 = x.std::array<T, N>::begin();

  while (first != last) {
    (*first++) += d * (*first2++) + e * (*first3++) + f * (*first4++);
  }
}

/**************************************************/

template<size_t N, typename T>
inline void
numfixarray<N, T>::sadd(double a, double d, const numfixarray<N, T>& v)
{
  nvp first = std::array<T, N>::begin();
  cnvp last = std::array<T, N>::end();
  cnvp first2 = v.std::array<T, N>::begin();

  while (first != last) {
    (*first) = a * (*first) + d * (*first2++);
    first++;
  }
}
} // namespace Gascoigne

#endif
