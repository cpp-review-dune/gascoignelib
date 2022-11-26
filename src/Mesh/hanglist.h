/**
 *
 * Copyright (C) 2004, 2005, 2008, 2011 by the Gascoigne 3D authors
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

#ifndef __hanglist_h
#define __hanglist_h

#include <string>

#include "../Interface/gascoignehash.h"

#include "edgearray.h"
#include "hang.h"

/*------------------------------------------------------*/

namespace Gascoigne {

//
/// This hash function has to be consistent with the operator "=="
/// for EdgeArrays, i.e. permutated fixarrays
//

template<int N>
class EdgeHash
{
public:
  size_t operator()(const EdgeArray<N>& h) const
  {
    size_t hash = std::hash<IndexType>{}(h[0]);
    for (size_t i = 1; i < N; ++i) {
      hash ^= std::hash<IndexType>{}(h[i]) << i;
    }
    return hash;
  }
};

/*------------------------------------------------------*/

template<int N>
class HangList : public std::unordered_map<EdgeArray<N>, Hang, EdgeHash<N>>
{
protected:
public:
  typedef typename std::unordered_map<EdgeArray<N>, Hang, EdgeHash<N>>::iterator
    iterator;
  typedef
    typename std::unordered_map<EdgeArray<N>, Hang, EdgeHash<N>>::const_iterator
      const_iterator;

  void update(const std::vector<IndexType>&);
  void update(const std::vector<IndexType>&, const std::vector<IndexType>&);
  void make_consistent(HangList<N>&);
  void move(HangList<N>& src, iterator& p);
  HangList<N>& operator=(const HangList<N>& A);
  void BinWrite(std::ostream& s) const;
  void BinRead(std::istream& s);
};

/*------------------------------------------------------*/

template<int N>
std::ostream&
operator<<(std::ostream& s, const HangList<N>& A);

template<int N>
std::istream&
operator>>(std::istream& s, HangList<N>& A);
} // namespace Gascoigne

/*------------------------------------------------------*/

#endif
