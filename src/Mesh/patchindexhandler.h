/**
 *
 * Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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

#ifndef __PatchIndexHandler_h
#define __PatchIndexHandler_h

#include "../Interface/gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne {
class PatchIndexHandler
{
protected:
  bool haspatch, hasq4patch;
  nvector<IndexVector> indexofpatch, indexofq4patch;
  nvector<IndexVector> patch2cell, q4patch2cell;

  int dim;

public:
  PatchIndexHandler()
    : haspatch(false)
    , hasq4patch(false)
  {}
  ~PatchIndexHandler() {}

  int& GetDim() { return dim; }
  bool& GetHasPatch() { return haspatch; }
  bool& GetHasQ4Patch() { return hasq4patch; }
  nvector<IndexVector>& GetIndex() { return indexofpatch; }
  const nvector<IndexVector>& GetIndex() const { return indexofpatch; }
  nvector<IndexVector>& GetIndexQ4() { return indexofq4patch; }
  const nvector<IndexVector>& GetIndexQ4() const
  {
    assert(hasq4patch);
    return indexofq4patch;
  }

  IndexVector& GetPatch2Cell(size_t i)
  {
    assert(i < patch2cell.size());
    return patch2cell[i];
  }
  const IndexVector& GetPatch2Cell(size_t i) const
  {
    assert(i < patch2cell.size());
    return patch2cell[i];
  }
  IndexVector& GetQ4Patch2Cell(size_t i)
  {
    assert(i < q4patch2cell.size());
    return q4patch2cell[i];
  }
  const IndexVector& GetQ4Patch2Cell(size_t i) const
  {
    assert(hasq4patch && i < q4patch2cell.size());
    return q4patch2cell[i];
  }

  nvector<IndexVector>& GetAllPatch2Cell() { return patch2cell; }
  const nvector<IndexVector>& GetAllPatch2Cell() const { return patch2cell; }
  nvector<IndexVector>& GetAllQ4Patch2Cell() { return q4patch2cell; }
  const nvector<IndexVector>& GetAllQ4Patch2Cell() const
  {
    assert(hasq4patch);
    return q4patch2cell;
  }

  size_t npatches() const { return indexofpatch.size(); }
  size_t nq4patches() const { return indexofq4patch.size(); }
  bool HasPatch() const { return haspatch; }
  bool HasQ4Patch() const { return hasq4patch; }
  int Dim() const { return dim; }

  const IndexVector& IndicesOfPatch(int i) const { return indexofpatch[i]; }
  const IndexVector& IndicesOfQ4Patch(int i) const
  {
    assert(hasq4patch);
    return indexofq4patch[i];
  }
  IndexVector Q2IndicesOfQ4Patch(int i) const;
  IndexVector CoarseIndices(int iq) const;
  IndexVector CoarseIndicesQ4(int iq) const;

  int nodes_per_patch() const
  {
    if (dim == 2)
      return 9;
    return 27;
  }
  int nodes_per_q4patch() const
  {
    if (dim == 2)
      return 25;
    return 125;
  }
};
} // namespace Gascoigne

#endif
