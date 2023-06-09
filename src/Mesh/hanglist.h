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

#include "edgearray.h"
#include "hang.h"
#include <string>

#include "gascoignehash.h"

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
  int operator()(const EdgeArray<N>& h) const { return h.sum(); }
};

#define HANGMAP HASHMAP<EdgeArray<N>, Hang, EdgeHash<N>>

/*------------------------------------------------------*/

template<int N>
class HangList : public HANGMAP
{
protected:
public:
  typedef typename HANGMAP::iterator iterator;
  typedef typename HANGMAP::const_iterator const_iterator;

  void update(const std::vector<int>&);
  void update(const std::vector<int>&, const std::vector<int>&);
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

#undef HANGMAP

#endif
