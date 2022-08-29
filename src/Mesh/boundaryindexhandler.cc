/**
 *
 * Copyright (C) 2004, 2005, 2007, 2009, 2011 by the Gascoigne 3D authors
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

#include "boundaryindexhandler.h"
#include "set2vec.h"
#include "stlio.h"
#include <set>

using namespace std;

/*--------------------------------------------------------------*/

namespace Gascoigne {
ostream&
operator<<(ostream& s, const BoundaryIndexHandler& A)
{
  const IndexSet& colors = A.GetColors();
  cerr << "All Colors: " << colors << endl;
  for (auto p = colors.begin(); p != colors.end(); ++p) {
    cerr << "color: " << *p;
    cerr << "\n\tVertices: " << A.Verteces(*p);
    cerr << "\n\tCells: " << A.Cells(*p);
    cerr << "\n\tLocalind: " << A.Localind(*p);
    cerr << endl;
  }
  return s;
}

/*--------------------------------------------------------------*/

void
BoundaryIndexHandler::check() const
{
  const IndexSet& colors = GetColors();
  cerr << "All Colors: " << colors << endl;
  for (auto p = colors.begin(); p != colors.end(); ++p) {
    cerr << "color: " << *p;
    const IndexVector& v = Verteces(*p);
    for (IndexType i = 0; i < v.size(); i++) {
      if (v[i] < 0) {
        cerr << "....BoundaryIndexHandler::check() ERROR\n";
        abort();
      }
    }
  }
  cerr << endl;
}

/*--------------------------------------------------------------*/

void
BoundaryIndexHandler::Equal(const IndexSet& col,
                            const VecMap& v,
                            const VecMap& c,
                            const VecMap& l,
                            const VecMap& patch,
                            const VecMap& patchlocal)
{
  AllColors = col;
  verteces = v;
  cells = c;
  localci = l;
  patches = patch;
  localpi = patchlocal;
}

/*--------------------------------------------------------------*/

void
BoundaryIndexHandler::CopySetToVector(const vector<IndexSet>& H,
                                      const IndexVector& colorvec,
                                      VecMap& dst) const
{
  for (IndexType i = 0; i < H.size(); i++) {
    IndexVector v;
    Set2Vec(v, H[i]);
    IndexType color = colorvec[i];
    dst.insert(make_pair(color, v));
  }
}

/*--------------------------------------------------------------*/

void
BoundaryIndexHandler::clear()
{
  AllColors.clear();
  verteces.clear();
  cells.clear();
  localci.clear();
}

/*--------------------------------------------------------------*/

const IndexVector&
BoundaryIndexHandler::Verteces(IndexType color) const
{
  VecMap::const_iterator p = verteces.find(color);
  if (p == verteces.end()) {
    cerr << "BoundaryIndexHandler::Vertices\tcolor not found: " << color
         << endl;
    abort();
  }
  return p->second;
}

/*--------------------------------------------------------------*/

const IndexVector&
BoundaryIndexHandler::Cells(IndexType color) const
{
  VecMap::const_iterator p = cells.find(color);
  if (p == cells.end()) {
    cerr << "BoundaryIndexHandler::Cells\tcolor not found: " << color << endl;
    abort();
  }
  return p->second;
}

/*--------------------------------------------------------------*/

const IndexVector&
BoundaryIndexHandler::Localind(IndexType color) const
{
  VecMap::const_iterator p = localci.find(color);
  if (p == localci.end()) {
    cerr << "BoundaryIndexHandler::Localind\tcolor not found: " << color
         << endl;
    abort();
  }
  return p->second;
}

/*--------------------------------------------------------------*/

const IndexVector&
BoundaryIndexHandler::Patches(IndexType color) const
{
  VecMap::const_iterator p = patches.find(color);
  if (p == patches.end()) {
    cerr << "BoundaryIndexHandler::Patches\tcolor not found: " << color << endl;
    abort();
  }
  return p->second;
}

/*--------------------------------------------------------------*/

const IndexVector&
BoundaryIndexHandler::LocalPatchind(IndexType color) const
{
  VecMap::const_iterator p = localpi.find(color);
  if (p == localpi.end()) {
    cerr << "BoundaryIndexHandler::LocalPatchind\tcolor not found: " << color
         << endl;
    abort();
  }
  return p->second;
}

/*--------------------------------------------------------------*/

void
BoundaryIndexHandler::SetPeriodicPairs(
  std::map<IndexType, std::map<IndexType, IndexType>> mm_PeriodicPairs)
{
  _PeriodicPairs = mm_PeriodicPairs;
}

/*--------------------------------------------------------------*/

const std::map<IndexType, std::map<IndexType, IndexType>>
BoundaryIndexHandler::GetPeriodicPairs() const
{
  return _PeriodicPairs;
}

} // namespace Gascoigne
