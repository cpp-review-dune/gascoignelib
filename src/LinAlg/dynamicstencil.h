/**
 *
 * Copyright (C) 2008 by the Gascoigne 3D authors
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

#ifndef __dynamicstencil_h
#define __dynamicstencil_h

#include <cassert>
#include <list>
#include <vector>

#include "../Interface/gascoigne.h"
#include "../Interface/stencilinterface.h"

namespace Gascoigne {

class DynamicStencil : public StencilInterface
{
protected:
  typedef std::list<IndexType>::const_iterator const_citerator;
  typedef std::list<IndexType>::iterator citerator;

  // the column for every entry in the row

public:
  std::vector<std::list<IndexType>> cols;

  IndexType n() const { return cols.size(); }

  const_citerator cstart(IndexType i) const
  {
    assert(i < n());
    return cols[i].begin();
  }
  const_citerator cstop(IndexType i) const
  {
    assert(i < n());
    return cols[i].end();
  }
  citerator cstart(IndexType i)
  {
    assert(i < n());
    return cols[i].begin();
  }
  citerator cstop(IndexType i)
  {
    assert(i < n());
    return cols[i].end();
  }
};
} // namespace Gascoigne

#endif
