/**
 *
 * Copyright (C) 2004, 2005, 2008 by the Gascoigne 3D authors
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

#ifndef __sparsestructure_h
#define __sparsestructure_h

#include <vector>

#include "../Interface/gascoigne.h"
#include "../Interface/sparsestructureinterface.h"

/*------------------------------------------------------------------------*/

namespace Gascoigne {
class SparseStructure : public SparseStructureInterface
{
protected:
  typedef std::vector<IndexSet> Indices;

  IndexType sntot;
  Indices sindices;

public:
  typedef IndexSet::iterator iterator;
  typedef IndexSet::const_iterator const_iterator;

  SparseStructure()
    : sntot(0)
    , sindices(0)
  {}

  void clear()
  {
    sntot = 0;
    sindices.clear();
  }

  IndexType n() const { return sindices.size(); }
  IndexType ntotal() const { return sntot; }

  const Indices& indices() const { return sindices; }
  const IndexSet& row(IndexType i) const
  {
    assert((i >= 0) && (i < sindices.size()));
    return sindices[i];
  }
  IndexSet& row(IndexType i)
  {
    assert((i >= 0) && (i < sindices.size()));
    return sindices[i];
  }
  IndexSet::iterator rowbegin(IndexType i)
  {
    assert((i >= 0) && (i < sindices.size()));
    return row(i).begin();
  }
  IndexSet::iterator rowend(IndexType i)
  {
    assert((i >= 0) && (i < sindices.size()));
    return row(i).end();
  }
  IndexSet::const_iterator rowbegin(IndexType i) const
  {
    assert((i >= 0) && (i < sindices.size()));
    return row(i).begin();
  }
  IndexSet::const_iterator rowend(IndexType i) const
  {
    assert((i >= 0) && (i < sindices.size()));
    return row(i).end();
  }
  IndexType rowsize(IndexType i) const
  {
    assert((i >= 0) && (i < sindices.size()));
    return row(i).size();
  }

  friend std::ostream& operator<<(std::ostream& s, const SparseStructure& A);
  void statistics(std::ostream&) const;

  SparseStructure& operator=(const SparseStructure& B);
  void build_begin(IndexType n);
  void build_clear(IndexType i);
  void build_add(IndexType i, IndexType j) { row(i).insert(j); }
  template<class IT>
  void build_add(IndexType i, IT lsf, IT lsl)
  {
    for (IT p = lsf; p != lsl; p++)
      row(i).insert(*p);
  }
  template<class IT>
  void build_add(IT lsf, IT lsl)
  {
    for (IT p = lsf; p != lsl; p++)
      build_add(*p, lsf, lsl);
  }
  template<class IT>
  void build_add(IT rf, IT rl, IT cf, IT cl)
  {
    for (IT p = rf; p != rl; p++)
      build_add(*p, cf, cl);
  }
  void build_end();

  void hanging_node(IndexType, IndexType, IndexType);

  void enlarge(const SparseStructure&);
  void enlarge_lu();
  void enlarge_for_lu(const IntVector& perm);
};
} // namespace Gascoigne

#endif
