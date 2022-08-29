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

#ifndef __vecalgo_h
#define __vecalgo_h

#include <set>
#include <vector>

#include "../Interface/gascoigne.h"

namespace Gascoigne {
void
transfer(IndexType n,
         std::vector<IndexType>& tr,
         const std::set<IndexType>& del);
void
transfer(IndexType n, std::vector<IndexType>& tr, std::vector<IndexType>& del);

template<class C>
void
compress(std::vector<C>& dst, const std::vector<IndexType>& src)
{
  // IndexType n = 0;
  IndexType mmax = 0;

  for (IndexType i = 0; i < src.size(); i++) {
    IndexType j = src[i];
    if (j >= 0) {
      dst[j] = dst[i];
      mmax = std::max(mmax, j);
      // n++;
    }
  }
  // dst.resize(n);
  dst.resize(mmax + 1);
}
} // namespace Gascoigne

#endif
