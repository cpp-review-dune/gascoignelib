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

#ifndef __edgemanager_h
#define __edgemanager_h

#include "../Interface/gascoigne.h"

#include "edge.h"
#include "hangcontainer2d.h"
#include "quadlawandorder.h"

/*---------------------------------------------------*/

namespace Gascoigne {
class EdgeManager
{

protected:
  std::vector<Edge>& edges;
  std::vector<Quad>& quads;
  const IndexVector& co2n;
  IndexVector& eo2n;

  IndexVector SwappedEdge;
  QuadLawAndOrder QuadLaO;

  void Update();
  void InnerEdges(const IndexSet& CellRefList);
  void OuterEdges(const HangContainer2d& hangset);
  void OldHangings(HangContainer2d& hangset, const IndexSet& CellRefList);
  void SwappedEdges();
  void NeighbourTester() const;

  void BSETest() const;

public:
  EdgeManager(std::vector<Edge>&,
              std::vector<Quad>&,
              const IndexVector& con,
              IndexVector& eon);

  const Quad& quad(IndexType i) const { return quads[i]; }
  Quad& quad(IndexType i) { return quads[i]; }

  std::array<IndexType, 2> ChildrenOfEdge(IndexType e) const;

  bool EdgeIsHanging(IndexType e) const;
  bool EdgeIsHanging(const Edge& e) const;

  void LoadEdgeElimination(IndexVector& edel,
                           const IndexSet& CellCoarseList,
                           const HangContainer2d& hangset) const;
  void Build(const IndexSet& CellRefList, HangContainer2d&);
  void DeleteEdges();
  void InitEdges();
  void SortHangings();
};
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif
