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


#ifndef __fixarray_h
#define __fixarray_h

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <stdlib.h>

/*-------------------------------------------------*/

namespace Gascoigne
{
  template <int N, class T>
  class fixarray : public std::array<T, N>
  {

  public:
  protected:
    void array_copy(const_iterator q)
    {
      iterator p(begin());
      const_iterator pe(end());
      while (p != pe)
        *p++ = *q++;
    }

  public:
    std::array< T,N>()
        : std::array<T, N>()
    {
      BasicInit(T());
    }
    std::array< T,N>(const T &d)
    {
      BasicInit(d);
    }
    std::array< T,N>(const std::array< T,N> &v)
    {
      BasicInit(T());
      array_copy(v.begin());
      // copy(v.begin(),v.end(),begin());
    }
    fixarray(const_iterator b)
    {
      BasicInit(T());
      array_copy(b);
    }

    virtual ~fixarray()
    {
      //       Destroy(begin(),end());
    }

    void BasicInit(const T &d)
    {
      // Braucht man das wirklich ???
      //       for(int i=0;i<N;i++)  construct(&(val[i]),d);
      for (int i = 0; i < N; i++)
        val[i] = d;
    }

    const T *begin() const
    {
      return &(val[0]);
    }
    const T *end() const
    {
      return &(val[0]) + N;
    }
    T *begin()
    {
      return &(val[0]);
    }
    T *end()
    {
      return &(val[0]) + N;
    }

    size_t size() const
    {
      return N;
    }
    const T &operator[](int i) const
    {
      return val[i];
    }
    T &operator[](int i)
    {
      return val[i];
    }

    std::array< T,N> &operator=(const T &d)
    {
      iterator p(end());
      while (p > begin())
        *--p = d;
      return *this;
    }

    std::array< T,N> &operator=(const std::array< T,N> &v)
    {
      iterator p(begin());
      const_iterator q(v.begin());
      while (p < end())
        *p++ = *q++;
      return *this;
    }

    bool operator<(const std::array< T,N> &v) const
    {
      const_iterator p(begin());
      const_iterator q(v.begin());
      while (p < end())
      {
        if (*p < *q)
          return 1;
        if (*q < *p)
          return 0;
        p++;
        q++;
      }
      return 0;
    }
    bool operator!=(const std::array< T,N> &v) const
    {
      const_iterator p(begin());
      const_iterator q(v.begin());
      while (p < end())
      {
        if (*p != *q)
          return 1;
        p++;
        q++;
      }
      return 0;
    }


    std::ostream &put(std::ostream &s) const
    {
      copy(begin(), end(), std::ostream_iterator<T>(s, " "));
      return s;
    }

    std::istream &get(std::istream &s)
    {
      typename std::array< T,N>::iterator p;
      for (p = begin(); p != end(); p++)
        s >> *p;
      return s;
    }

    void read_data(std::istream &s)
    {
      size_t n;
      s >> n;
      if (size() != n)
      {
        std::cerr << "read_data(): wrong size in fixarray" << N << " " << n
                  << std::endl;
        exit(1);
      }
      s >> *this;
    }

    void write_data(std::ostream &s) const
    {
      s << size() << std::endl;
      s << *this;
    }

    void BinWrite(std::ostream &s) const
    {
      int sizeT = sizeof(T);
      for (int i = 0; i < N; i++)
      {
        s.write(reinterpret_cast<const char *>(&(operator[](i))), sizeT);
      }
    }

    void BinRead(std::istream &s)
    {
      int sizeT = sizeof(T);
      for (int i = 0; i < N; i++)
      {
        s.read(reinterpret_cast<char *>(&(operator[](i))), sizeT);
      }
    }
  };

  /*-------------------------------------------------*/

  template <int N, class T>
  bool operator==(const std::array< T,N> &x, const std::array< T,N> &y)
  {
    return std::equal(x.begin(), x.end(), y.begin());
  }

  /*-------------------------------------------------*/

  class fixarrayHash
  {
  public:
    template <int N, class T>
    int operator()(const std::array< T,N> &h) const
    {
      return static_cast<int>(h[0]);
    }
  };

  template <int N, class T>
  std::ostream &operator<<(std::ostream &s, const std::array< T,N> &A);
  template <int N, class T>
  std::istream &operator>>(std::istream &s, std::array< T,N> &A);
}


#endif
