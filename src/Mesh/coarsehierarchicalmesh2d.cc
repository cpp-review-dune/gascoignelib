/**
 *
 * Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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

#include "coarsehierarchicalmesh2d.h"
#include "set2vec.h"

/*---------------------------------------------------*/

namespace Gascoigne {
CoarseHierarchicalMesh2d::CoarseHierarchicalMesh2d(const HierarchicalMesh2d& HM)
  : HierarchicalMesh2d(HM)
{}

/*---------------------------------------------------*/

void
CoarseHierarchicalMesh2d::BasicInit(int depth)
{
  assert(1 <= depth && depth <= 2);

  if (depth == 1) {
    loop(cn2o);
  } else if (depth == 2) {
    IndexVector cn2oA, cn2oB;
    IndexVector co2n2;
    loop(cn2oA);
    co2n2 = co2n;
    loop(cn2oB);
    IndexVector co2n3(co2n2.size());
    for (size_t i = 0; i < co2n2.size(); i++) {
      if (co2n2[i] < 0) {
        co2n3[i] = co2n2[i];
      } else {
        assert(co2n2[i] < co2n.size());
        co2n3[i] = co2n[co2n2[i]];
      }
    }
    co2n.resize(co2n3.size());
    co2n = co2n3;

    cn2o.resize(ncells());
    cn2o = -1;
    for (size_t i = 0; i < cn2oB.size(); i++) {
      int j = cn2oB[i];
      if (j >= 0) {
        int k = cn2oA[j];
        if (k >= 0) {
          cn2o[i] = k;
        }
      }
    }
  }
}

/*---------------------------------------------------*/

void
CoarseHierarchicalMesh2d::loop(IndexVector& dst)
{
  global_coarse();
  dst.resize(ncells());
  dst = -1;
  for (size_t i = 0; i < co2n.size(); i++) {
    int j = co2n[i];
    assert(j < ncells());
    if (j >= 0)
      dst[j] = i;
  }
}

/*---------------------------------------------------*/

void
CoarseHierarchicalMesh2d::refine(const IndexVector& cell_ref_old,
                                 const IndexVector& cell_coarse_old)
{
  CellRefList.clear();
  CellCoarseList.clear();

  IndexVector cell_ref(0), cell_coarse(0);

  for (size_t i = 0; i < cell_ref_old.size(); i++) {
    int newc = co2n[cell_ref_old[i]];
    if (newc >= 0)
      cell_ref.push_back(newc);
  }
  for (size_t i = 0; i < cell_coarse_old.size(); i++) {
    int newc = co2n[cell_coarse_old[i]];
    if (newc >= 0)
      cell_coarse.push_back(newc);
  }
  _refine2d(CellRefList, CellCoarseList, cell_ref, cell_coarse);
}

/*---------------------------------------------------*/

void
CoarseHierarchicalMesh2d::GetRefinedList(IndexVector& ref)
{
  ref.resize(0);
  IndexVector ref2;
  Set2Vec(ref2, CellRefList);
  for (size_t i = 0; i < ref2.size(); i++) {
    int j = cn2o[ref2[i]];
    assert(j >= 0);
    ref.push_back(j);
  }
}

/*---------------------------------------------------*/

void
CoarseHierarchicalMesh2d::GetCoarsedList(IndexVector& coarse)
{
  coarse.resize(0);
  IndexVector coarse2;
  Set2Vec(coarse2, CellCoarseList);
  for (size_t i = 0; i < coarse2.size(); i++) {
    int j = cn2o[coarse2[i]];
    if (j >= 0)
      coarse.push_back(j);
  }
}
} // namespace Gascoigne

/*---------------------------------------------------*/
