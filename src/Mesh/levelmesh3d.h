/**
 *
 * Copyright (C) 2004, 2006, 2007 by the Gascoigne 3D authors
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

#ifndef __levelmesh3d_h
#define __levelmesh3d_h

#include "boundaryindexhandler.h"
#include "gascoigne.h"
#include "hierarchicalmesh3d.h"
#include "index.h"
#include "patchindexhandler.h"

/*--------------------------------------------------------------*/

namespace Gascoigne {
class LevelMesh3d : public Index
{
protected:
  typedef std::map<IndexType, std::array<IndexType, 3>> QuadraticHNStructure3;
  typedef std::map<IndexType, std::array<IndexType, 9>> QuadraticHNStructure9;

  typedef std::map<IndexType, std::array<IndexType, 6>> QuarticHNStructure5;
  typedef std::map<IndexType, std::array<IndexType, 26>> QuarticHNStructure25;

  const HierarchicalMesh3d* HMP;

  void check_leveljump() const;
  void fill_opis(IndexSet& dst, IndexSet& oldquads) const;
  void fill_enkel(IndexSet& dst, const Hex& Q) const;
  void fill_childs(IndexSet& dst, const Hex& Q) const;
  bool EnkelUniform(const Hex& Q) const;
  bool BuildFathers(IndexSet& Vaeter) const;
  IndexType hexindex(IndexType i) const { return Index::Hexl2g(i); }
  void InitCells(IndexType n);
  void InitNodes(IndexType n);
  void InitEdges(IndexType n);

  IndexType refine_level(IndexType n) const;
  void ConstructNodesOnFaceQ4(std::array<IndexType, 81>& nodesonface,
                              IndexType vater,
                              IndexType ni) const;
  void InsertHangingFacesQ4(QuarticHNStructure25& hnq4face,
                            const std::array<IndexType, 81>& nodesonface) const;
  void InsertHangingEdgesQ4(QuarticHNStructure5& hnq4,
                            const std::array<IndexType, 81>& nodesonface) const;
  void ConstructNodesOnFace(std::array<IndexType, 25>& nodesonface,
                            IndexType vater,
                            IndexType ni) const;
  void InsertHangingFaceQ4(QuarticHNStructure25& hnq4face,
                           const std::array<IndexType, 81>& nodesonface,
                           IndexType n1,
                           IndexType n2,
                           IndexType n3,
                           IndexType n4,
                           const std::array<IndexType, 25>& I) const;
  void InsertHangingEdgeQ4(QuarticHNStructure5& hnq4,
                           const std::array<IndexType, 81>& nof,
                           IndexType n1,
                           IndexType n2,
                           IndexType n3,
                           IndexType n4,
                           IndexType i1,
                           IndexType i2,
                           IndexType i3,
                           IndexType i4,
                           IndexType i5) const;

public:
  LevelMesh3d(const HierarchicalMesh* hmp);
  ~LevelMesh3d();

  IndexType ncells() const { return Index::HexSize(); }
  IndexType dimension() const { return HMP->dimension(); }

  IndexType vertex_of_cell(IndexType i, IndexType ii) const
  {
    return Index::Vertexg2l(HMP->vertex_of_cell(hexindex(i), ii));
  }

  const Vertex3d& vertex3d(IndexType i) const
  {
    return HMP->vertex3d(Index::Vertexl2g(i));
  }

  /*----- Functions -----*/

  const Hex& hex(IndexType i) const { return HMP->hex(hexindex(i)); }

  void BasicInit(const IndexSet& n, const IndexSet& o);

  /*----- Functions fuer Patch-----*/

  void construct_lists(IndexSet& newquads, IndexSet& oldquads) const;

  void ConstructHangingStructureQuadratic(
    QuadraticHNStructure3& hnq2,
    QuadraticHNStructure9& hnq2face) const;
  void ConstructHangingStructureQuartic(QuarticHNStructure5& hnq4,
                                        QuarticHNStructure25& hnq4face) const;

  void InitBoundaryHandler(BoundaryIndexHandler& BI,
                           const PatchIndexHandler& PIH) const;
  bool ConstructCellIndOfPatch(IndexVector& dstc) const;
  void ConstructIndOfPatch(nvector<IndexVector>& dstv) const;
};
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif
