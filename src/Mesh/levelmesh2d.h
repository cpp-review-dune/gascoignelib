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

#ifndef __levelmesh2d_h
#define __levelmesh2d_h

#include "boundaryindexhandler.h"
#include "hierarchicalmesh2d.h"
#include "index.h"
#include "patchindexhandler.h"

/*--------------------------------------------------------------*/

namespace Gascoigne {
class LevelMesh2d : public Index
{
protected:
  typedef std::map<IndexType, std::array<IndexType, 3>> QuadraticHNStructure3;
  typedef std::map<IndexType, std::array<IndexType, 6>> QuarticHNStructure5;

  const HierarchicalMesh2d* HMP;

  void check_leveljump() const;
  void fill_opis(IndexSet& dst, IndexSet& oldquads) const;
  void fill_enkel(IndexSet& dst, const Quad& Q) const;
  void fill_childs(IndexSet& dst, const Quad& Q) const;
  bool EnkelUniform(const Quad& Q) const;
  bool BuildFathers(std::set<IndexType>& Vaeter) const;
  void InitCells(IndexType n);
  void InitNodes(IndexType n);
  void InitEdges(IndexType n);

public:
  /*----- Constructor -----*/

  LevelMesh2d(const HierarchicalMesh* hmp);
  ~LevelMesh2d();

  const HierarchicalMesh2d* GetHierarchicalMesh() const { return HMP; }

  IndexType ncells() const { return Index::QuadSize(); }
  const Quad& quad(IndexType i) const { return HMP->quad(Quadl2g(i)); }
  const Vertex2d& vertex2d(IndexType i) const
  {
    return HMP->vertex2d(Vertexl2g(i));
  }

  IndexType vertex_of_cell(IndexType i, IndexType j) const
  {
    return Vertexg2l(HMP->vertex_of_cell(Quadl2g(i), j));
  }

  /*----- Functions -----*/

  bool EdgeIsHangingGlobalIndex(IndexType i) const;

  void BasicInit(const IndexSet& n, const IndexSet& o);

  /*----- Functions for patches -----*/

  void construct_lists(IndexSet& newquads, IndexSet& oldquads) const;
  void ConstructHangingStructureQuadratic(QuadraticHNStructure3& hnq2) const;
  void ConstructHangingStructureQuartic(QuarticHNStructure5& hnq4) const;
  void InitBoundaryHandler(BoundaryIndexHandler& BI,
                           const PatchIndexHandler& PIH) const;
  void ConstructIndOfPatch(nvector<IndexVector>& dstv) const;
  bool ConstructCellIndOfPatch(IndexVector& dstc) const;
};
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif
