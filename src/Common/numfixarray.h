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
#include <cstring>
#include <array>
#include "nvector.h"

/*-------------------------------------------------*/

namespace Gascoigne
{
template <size_t N, typename T>
class numfixarray : public std::array<T, N>
{
  typedef typename std::array<T, N>::iterator nvp;
  typedef typename std::array<T, N>::const_iterator cnvp;
  static constexpr auto ALIGN = N * sizeof(double);

public:
  ~numfixarray()
  {
  }
  numfixarray()
  {
    std::array<T, N>::fill(0.0);
  }
  /*   numfixarray(size_t n)               : std::array<T,N>(n)   {} */
  numfixarray(const T& d)
  {
    std::array<T, N>::fill(d);
  }
  numfixarray(const numfixarray<N, T>& v) : std::array<T, N>(v)
  {
    std::memcpy(this->data(), v.data(), N * sizeof(double));
  }

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

  int n() const
  {
    return std::array<T, N>::size();
  }
  void zero();
  void equ(const T&);
  void equ(double, const numfixarray&);
  void equ(double, const numfixarray&, double, const numfixarray&);
  void equ(double, const numfixarray&, double, const numfixarray&, double,
           const numfixarray&);
  void equ(double, const numfixarray&, double, const numfixarray&, double,
           const numfixarray&, double, const numfixarray&);
  void equ(double, const numfixarray&, double, const numfixarray&, double,
           const numfixarray&, double, const numfixarray&, double, const numfixarray&,
           double, const numfixarray&);
  void sequ(double, double, const numfixarray&);
  void add(double);
  void add(T, const numfixarray&);
  void add(double, const numfixarray&, double, const numfixarray&);
  void add(double, const numfixarray&, double, const numfixarray&, double,
           const numfixarray&);
  void sadd(double, double, const numfixarray&);
  double max() const;
  double norm() const;
  double norm_l8() const;
  double norm_l1() const;
  double norm_l2() const
  {
    return norm();
  }

  void normalise()
  {
    T d = 0.;
    for (int i = 0; i < n(); i++)
    {
      d += (*this)[i] * (*this)[i];
    }
    d = 1. / sqrt(d);
    for (int i = 0; i < n(); i++)
    {
      (*this)[i] *= d;
    }
  }
};

/**************************************************/

template <size_t N, typename T>
inline double numfixarray<N, T>::norm_l8() const
{
  cnvp first = std::array<T, N>::cbegin();
  cnvp last  = std::array<T, N>::cend();

  double d = 0.;
  while (first != last)
  {
    d = std::max(d, fabs(*first));
    first++;
  }
  return d;
}

/**************************************************/

template <size_t N, typename T>
inline double numfixarray<N, T>::norm() const
{
  // cnvp first = std::array<T, N>::cbegin();
  // cnvp last  = std::array<T, N>::cend();

  double n = 0.;
  // while (first != last)
  // {
  //   n += ((*first)) * ((*first));
  //   first++;
  // }
  const double* __restrict__ tc = static_cast<const double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    n += tc[i] * tc[i];
  }
  return sqrt(n);
}

/**************************************************/

template <size_t N, typename T>
inline double numfixarray<N, T>::norm_l1() const
{
  // cnvp first = std::array<T, N>::cbegin();
  // cnvp last  = std::array<T, N>::cend();

  double n = 0.;
  // while (first != last)
  // {
  //   n += fabs((*first++));
  // }
  const double* __restrict__ tc = static_cast<const double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    n += fabs(tc[i]);
  }
  return n;
}

/**************************************************/

template <size_t N, typename T>
inline numfixarray<N, T>& numfixarray<N, T>::operator=(const T& d)
{
  // nvp first = std::array<T, N>::begin();
  // cnvp last = std::array<T, N>::cend();

  // while (first != last)
  // {
  //   *first++ = d;
  // }
  double* __restrict__ tc = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = d;
  }
  return *this;
}

/**************************************************/

template <size_t N, typename T>
inline numfixarray<N, T>& numfixarray<N, T>::operator=(const numfixarray<N, T>& v)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp vfirst = v.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   *first++ = *vfirst++;
  // }
  // return *this;
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = vc[i];
  }
  return *this;
}

/**************************************************/

template <size_t N, typename T>
inline numfixarray<N, T>& numfixarray<N, T>::operator*=(const numfixarray<N, T>& d)
{
  // nvp first = std::array<T, N>::begin();
  // cnvp last = std::array<T, N>::cend();
  // cnvp fd   = d.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first++) *= (*fd++);
  // }
  // return *this;
  std::array<double, N>* __restrict__ dc = static_cast<std::array<double, N>*>(&d);
  std::array<double, N>* __restrict__ tc = static_cast<std::array<double, N>*>(&this);
  for (int i = 0; i < N; ++i)
  {
    (*tc)[i] *= (*dc)[i];
  }
  return *this;
}

/**************************************************/

template <size_t N, typename T>
inline numfixarray<N, T>& numfixarray<N, T>::operator*=(double d)
{
  // nvp first = std::array<T, N>::begin();
  // cnvp last = std::array<T, N>::cend();

  // while (first != last)
  // {
  //   (*first++) *= d;
  // }
  // return *this;
  double* __restrict__ tc = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] *= d;
  }
  return *this;
}

/**************************************************/

template <size_t N, typename T>
inline numfixarray<N, T>& numfixarray<N, T>::operator/=(double d)
{
  // nvp first = std::array<T, N>::begin();
  // cnvp last = std::array<T, N>::cend();

  // while (first != last)
  // {
  //   (*first++) /= d;
  // }
  double* __restrict__ tc = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] /= d;
  }
  return *this;
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::zero()
{
  // for (auto& el : *this)
  // {
  //   el = 0.;
  // }
  double* __restrict__ tc = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = 0.0;
  }
}

/**************************************************/

template <size_t N, typename T>
inline double numfixarray<N, T>::max() const
{
  double d   = -1.;
  cnvp first = std::array<T, N>::cbegin();
  cnvp last  = std::array<T, N>::cend();

  while (first != last)
  {
    d = std::max(d, fabs((*first)));
    first++;
  }
  return d;
}

/**************************************************/

template <size_t N, typename T>
inline double numfixarray<N, T>::operator*(const numfixarray<N, T>& v) const
{
  // cnvp first  = std::array<T, N>::cbegin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();

  double d = 0.;
  // while (first != last)
  // {
  //   d += (*first++) * (*first2++);
  // }

  const double* __restrict__ vc = static_cast<const double*>(v.data());
  const double* __restrict__ tc = static_cast<const double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    d += tc[i] * vc[i];
  }
  return d;
}

/**************************************************/

template <size_t N, typename T>
inline double numfixarray<N, T>::operator*(const nvector<T>& v) const
{
  // cnvp first  = std::array<T, N>::cbegin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.cbegin();

  double d = 0.;
  // while (first != last)
  // {
  //   d += (*first++) * (*first2++);
  // }
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    d += tc[i] * vc[i];
  }
  return d;
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::equ(const T& d)
{
  // nvp first = std::array<T, N>::begin();
  // cnvp last = std::array<T, N>::cend();

  // while (first != last)
  // {
  //   (*first++) = d;
  // }
  double* __restrict__ tc = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = d;
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::equ(double d, const numfixarray<N, T>& v)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first++) = d * (*first2++);
  // }
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = d * vc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::equ(double d, const numfixarray<N, T>& v, double e,
                                   const numfixarray<N, T>& w)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();
  // cnvp first3 = w.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first++) = d * (*first2++) + e * (*first3++);
  // }
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  const double* __restrict__ wc = static_cast<const double*>(w.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = d * vc[i] + e * wc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::equ(double d, const numfixarray<N, T>& v, double e,
                                   const numfixarray<N, T>& w, double f,
                                   const numfixarray<N, T>& x)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();
  // cnvp first3 = w.std::array<T, N>::cbegin();
  // cnvp first4 = x.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first++) = d * (*first2++) + e * (*first3++) + f * (*first4++);
  // }
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  const double* __restrict__ wc = static_cast<const double*>(w.data());
  const double* __restrict__ xc = static_cast<const double*>(x.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = d * vc[i] + e * wc[i] + f * xc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::equ(double d, const numfixarray<N, T>& v, double e,
                                   const numfixarray<N, T>& w, double f,
                                   const numfixarray<N, T>& x, double g,
                                   const numfixarray<N, T>& y)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();
  // cnvp first3 = w.std::array<T, N>::cbegin();
  // cnvp first4 = x.std::array<T, N>::cbegin();
  // cnvp first5 = y.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first++) = d * (*first2++) + e * (*first3++) + f * (*first4++) + g * (*first5++);
  // }
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  const double* __restrict__ wc = static_cast<const double*>(w.data());
  const double* __restrict__ xc = static_cast<const double*>(x.data());
  const double* __restrict__ yc = static_cast<const double*>(y.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = d * vc[i] + e * wc[i] + f * xc[i] + g * yc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::equ(double d, const numfixarray<N, T>& v, double e,
                                   const numfixarray<N, T>& w, double f,
                                   const numfixarray<N, T>& x, double g,
                                   const numfixarray<N, T>& y, double h,
                                   const numfixarray<N, T>& z, double i,
                                   const numfixarray<N, T>& zz)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();
  // cnvp first3 = w.std::array<T, N>::cbegin();
  // cnvp first4 = x.std::array<T, N>::cbegin();
  // cnvp first5 = y.std::array<T, N>::cbegin();
  // cnvp first6 = z.std::array<T, N>::cbegin();
  // cnvp first7 = zz.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first++) = d * (*first2++) + e * (*first3++) + f * (*first4++) + g * (*first5++)
  //                + h * (*first6++) + i * (*first7++);
  // }
  const double* __restrict__ vc  = static_cast<const double*>(v.data());
  const double* __restrict__ wc  = static_cast<const double*>(w.data());
  const double* __restrict__ xc  = static_cast<const double*>(x.data());
  const double* __restrict__ yc  = static_cast<const double*>(y.data());
  const double* __restrict__ zc  = static_cast<const double*>(z.data());
  const double* __restrict__ zzc = static_cast<const double*>(zz.data());
  double* __restrict__ tc        = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = d * vc[i] + e * wc[i] + f * xc[i] + g * yc[i] + h * zc[i] + i * zzc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::sequ(double s, double d, const numfixarray<N, T>& v)
{
  // equ(s, *this);
  // add(d, v);
  // return;

  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first) = s * (*first) + d * (*first2++);
  //   first++;
  // }
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] = s * tc[i] + d * vc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::add(double d)
{
  double* __restrict__ tc = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] += d;
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::add(T d, const numfixarray<N, T>& v)
{
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] += d * vc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::add(double d, const numfixarray<N, T>& v, double e,
                                   const numfixarray<N, T>& w)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();
  // cnvp first3 = w.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first++) += d * (*first2++) + e * (*first3++);
  // }
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  const double* __restrict__ wc = static_cast<const double*>(w.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] += d * vc[i] + e * wc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::add(double d, const numfixarray<N, T>& v, double e,
                                   const numfixarray<N, T>& w, double f,
                                   const numfixarray<N, T>& x)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();
  // cnvp first3 = w.std::array<T, N>::cbegin();
  // cnvp first4 = x.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first++) += d * (*first2++) + e * (*first3++) + f * (*first4++);
  // }

  const double* __restrict__ vc = static_cast<const double*>(v.data());
  const double* __restrict__ wc = static_cast<const double*>(w.data());
  const double* __restrict__ xc = static_cast<const double*>(x.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] += d * vc[i] + e * wc[i] + f * xc[i];
  }
}

/**************************************************/

template <size_t N, typename T>
inline void numfixarray<N, T>::sadd(double a, double d, const numfixarray<N, T>& v)
{
  // nvp first   = std::array<T, N>::begin();
  // cnvp last   = std::array<T, N>::cend();
  // cnvp first2 = v.std::array<T, N>::cbegin();

  // while (first != last)
  // {
  //   (*first) = a * (*first) + d * (*first2++);
  //   first++;
  // }
  const double* __restrict__ vc = static_cast<const double*>(v.data());
  double* __restrict__ tc       = static_cast<double*>(this->data());
  for (int i = 0; i < N; ++i)
  {
    tc[i] += a * tc[i] + d * vc[i];
  }
}
}  // namespace Gascoigne

#endif
