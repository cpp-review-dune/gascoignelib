/**
 *
 * Copyright (C) 2004 by the Gascoigne 3D authors
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

#include "../Interface/gascoigne.h"

#include "set2vec.h"

using namespace std;

namespace Gascoigne {

/*---------------------------------------------------*/
template<typename T>
void
Set2Vec(vector<T>& v, const set<T>& h)
{
  v.resize(h.size());
  size_t j = 0;
  for (auto p = h.begin(); p != h.end(); p++) {
    v[j++] = *p;
  }
}
template void
Set2Vec<IndexType>(vector<IndexType>& v, const set<IndexType>& h);

/*---------------------------------------------------*/

template<typename T>
void
Vec2Set(set<T>& h, const vector<T>& v)
{
  h.clear();
  for (size_t i = 0; i < v.size(); i++) {
    h.insert(v[i]);
  }
}

template void
Vec2Set<IndexType>(set<IndexType>& h, const vector<IndexType>& v);

} // namespace Gascoigne
