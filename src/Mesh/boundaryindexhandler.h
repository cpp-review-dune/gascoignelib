/**
 *
 * Copyright (C) 2004, 2007, 2009 by the Gascoigne 3D authors
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

#ifndef __boundaryindexhandler_h
#define __boundaryindexhandler_h

#include <map>
#include <set>

#include "../Interface/gascoigne.h"

/*--------------------------------------------------------------*/

namespace Gascoigne {
class BoundaryIndexHandler
{
protected:
  IndexSet AllColors;
  VecMap verteces, cells, localci, patches, localpi;

  std::map<IndexType, std::map<IndexType, IndexType>> _PeriodicPairs;

public:
  void CopySetToVector(const std::vector<IndexSet>&,
                       const IndexVector&,
                       VecMap&) const;

  void clear();

  const IndexSet& GetColors() const { return AllColors; }
  const VecMap& GetVertex() const { return verteces; }
  const VecMap& GetCell() const { return cells; }
  const VecMap& GetLocal() const { return localci; }
  const VecMap& GetPatch() const { return patches; }
  const VecMap& GetLocalPatch() const { return localpi; }

  IndexSet& GetColors() { return AllColors; }
  VecMap& GetVertex() { return verteces; }
  VecMap& GetCell() { return cells; }
  VecMap& GetLocal() { return localci; }
  VecMap& GetPatch() { return patches; }
  VecMap& GetLocalPatch() { return localpi; }

  const IndexVector& Verteces(IndexType col) const;
  const IndexVector& Cells(IndexType col) const;
  const IndexVector& Localind(IndexType col) const;
  const IndexVector& Patches(IndexType col) const;
  const IndexVector& LocalPatchind(IndexType col) const;

  void SetPeriodicPairs(
    std::map<IndexType, std::map<IndexType, IndexType>> mm_PeriodicPairs);
  const std::map<IndexType, std::map<IndexType, IndexType>> GetPeriodicPairs()
    const;

  void Equal(const IndexSet& col,
             const VecMap& v,
             const VecMap& c,
             const VecMap& l,
             const VecMap& patch,
             const VecMap& patchlocal);
  void check() const;
  friend std::ostream& operator<<(std::ostream& s,
                                  const BoundaryIndexHandler& A);
};
} // namespace Gascoigne

/*--------------------------------------------------------------*/

#endif
