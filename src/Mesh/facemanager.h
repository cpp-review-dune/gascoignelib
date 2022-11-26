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

#ifndef __facemanager_h
#define __facemanager_h

#include "../Interface/gascoigne.h"

#include "edge.h"
#include "hangcontainer3d.h"
#include "hexlawandorder.h"

/*---------------------------------------------------*/

namespace Gascoigne {
class FaceManager
{
protected:
  std::vector<Edge>& edges;
  std::vector<Hex>& hexs;
  const IndexVector& co2n;
  IndexVector& eo2n;

  IndexVector SwappedEdge;
  HexLawAndOrder HexLaO;

  void Update();
  void InnerFaces(const IndexSet& CellRefList);
  void OuterFaces(const HangContainer3d& hangset);
  void OldHangings(HangContainer3d& hangset3d, const IndexSet& CellRefList);
  void SwappedFaces();
  void NeighbourTester() const;
  void FillNeighbourFaces(const Hex& M, const Hex& S, const FaceVector& Face);

public:
  FaceManager(std::vector<Edge>&,
              std::vector<Hex>&,
              const IndexVector& con,
              IndexVector& eon);

  const Hex& hex(int i) const { return hexs[i]; }
  Hex& hex(int i) { return hexs[i]; }

  std::array<int, 2> ChildrenOfFace(int e) const;

  bool EdgeIsHanging(int e) const;
  bool EdgeIsHanging(const Edge& e) const;

  void LoadFaceElimination(IndexVector& edel,
                           const IndexSet& CellCoarseList,
                           const HangContainer3d& hangset) const;
  void Build(const IndexSet& CellRefList, HangContainer3d&);
  void DeleteFaces();
  void InitFaces();
  void SortHangings();
  void Check(const HangContainer3d& hangset) const;
};
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif
