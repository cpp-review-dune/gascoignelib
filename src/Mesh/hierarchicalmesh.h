/**
 *
 * Copyright (C) 2004, 2005, 2006, 2011 by the Gascoigne 3D authors
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

#ifndef __hierarchicalmesh_h
#define __hierarchicalmesh_h

#include <utility>

#include "compvector.h"
#include "curvedshapes.h"
#include "edge.h"
#include "gascoigne.h"
#include "hanglist.h"
#include "paramfile.h"
#include "triple.h"
#include "vertex.h"
#include <string>

/*---------------------------------------------------*/

namespace Gascoigne {
class HierarchicalMesh
{
protected:
  /*  typedef  */

  typedef std::pair<IndexType, IndexType> pint;
  typedef triple<IndexType, IndexType, IndexType> tint;
  typedef std::array<IndexType, 2> EdgeVector;
  typedef std::array<IndexType, 4> FaceVector;
  typedef std::vector<Edge> EdgeVec;
  typedef IndexSet::iterator IntSetIt;
  typedef IndexSet::const_iterator IntSetCIt;

  /*  Data  */

  int mnlevels;
  IndexType pdepth, etapatcher;
  IndexVector vo2n, eo2n, co2n;
  EdgeVec edges;

  mutable IndexType _i_showoutput;

  void update_edges(IndexVector&);
  virtual IndexType FindPatchDepth() const = 0;
  virtual void FillVertexLevels(IndexVector& dst) const = 0;
  virtual void RefineCoarseNodes(IndexSet& dst,
                                 const IndexVector& refnodes,
                                 const IndexVector& vertexlevel) const = 0;
  virtual void VertexToCells(IndexVector& dst,
                             const IndexSet& src,
                             const IndexVector& vertexlevel) const = 0;
  virtual void VertexToCellsCoarsening(
    IndexVector& dst,
    const IndexSet& src,
    const IndexVector& vertexlevel) const = 0;

public:
  virtual ~HierarchicalMesh();

  IndexType withfaces;

  HierarchicalMesh();
  HierarchicalMesh(const HierarchicalMesh&);
  HierarchicalMesh& operator=(const HierarchicalMesh&);

  /*  Zugriff  */

  virtual IndexType nnodes() const = 0;
  IndexType nlevels() const { return 1 + mnlevels; }
  IndexType nedges() const { return edges.size(); }

  const IndexVector* Vertexo2n() const { return &vo2n; }
  const IndexVector* Edgeo2n() const { return &eo2n; }
  const IndexVector* Cello2n() const { return &co2n; }

  IndexType Vertexo2n(IndexType i) const
  {
    assert(i < vo2n.size());
    return vo2n[i];
  }
  IndexType Edgeo2n(IndexType i) const
  {
    assert(i < eo2n.size());
    return eo2n[i];
  }
  IndexType Cello2n(IndexType i) const
  {
    assert(i < co2n.size());
    return co2n[i];
  }

  const Edge& edge(IndexType i) const
  {
    assert(i < edges.size());
    return edges[i];
  }
  const EdgeVec& edge() const { return edges; }

  virtual IndexType child(IndexType i, IndexType ii) const = 0;
  virtual IndexType nchilds(IndexType i) const = 0;

  virtual IndexType level(IndexType i) const = 0;
  virtual bool sleep(IndexType i) const = 0;

  virtual IndexType Vater(const IndexType i) const
  {
    std::cerr << "\"HierarchicalMesh::Vater\" not written!" << std::endl;
    abort();
  }
  virtual IndexVector Nachkommen(const IndexType i) const
  {
    std::cerr << "\"HierarchicalMesh::Nachkommen\" not written!" << std::endl;
    abort();
  }
  virtual IndexVector Geschwister(const IndexType i) const
  {
    std::cerr << "\"HierarchicalMesh::Geschwister\" not written!" << std::endl;
    abort();
  }
  virtual IndexVector Kinder(const IndexType i) const
  {
    std::cerr << "\"HierarchicalMesh::Kinder\" not written!" << std::endl;
    abort();
  }

  virtual IndexType nodes_per_cell(IndexType i) const = 0;
  virtual IndexType vertex_of_cell(IndexType i, IndexType ii) const = 0;

  void SetParameters(std::string gridname,
                     IndexType patchdepth,
                     IndexType epatcher);
  void SetParameters(IndexType patchdepth);
  void ReadFile(const std::string& gridname);
  void BasicInit(const ParamFile& pf, IndexType pdepth = 0);
  void global_refine(IndexType k);
  void global_patch_coarsen(IndexType k);
  void random_refine(double, IndexType k = 1);
  void random_patch_refine(double, IndexType k = 1);
  void random_patch_coarsen(double, IndexType k = 0);
  void random_double_patch_refine(double, IndexType k = 1);
  void clear_transfer_lists();
  virtual void write_gip(const std::string&) const = 0;
  virtual void write_gup(const std::string&) const = 0;
  virtual void write_inp(const std::string&) const = 0;

  virtual IndexType dimension() const { return 0; }
  virtual IndexType ncells() const = 0;
  IndexType patchdepth() const { return pdepth; }
  virtual IndexType nactivedescendants(IndexType i) const = 0;
  virtual IndexVector GetVertices(IndexType c) const = 0;

  bool CellIsCurved(IndexType iq) const
  {
    return GetBoundaryCellOfCurved(iq) != -1;
  }

  virtual std::set<IndexType> GetColors() const = 0;

  virtual void read_inp(const std::string&) = 0;
  virtual void read_gup(const std::string&) = 0;
  virtual void read_gip(const std::string&) = 0;
  virtual void refine(const IndexVector&, const IndexVector&) = 0;
  virtual void patch_refine(IndexVector&, IndexVector&) = 0;
  virtual void vertex_patch_refine(IndexVector& ref, IndexVector& coarse);
  virtual void vertex_patch_refine(IndexVector&);
  virtual void GetAwakePatchs(std::set<IndexType>&) const = 0;
  virtual void GetAwakeCells(std::set<IndexType>&) const = 0;
  virtual void ConstructQ2PatchMesh(IndexVector& pm) const = 0;
  virtual IndexVector ConstructQ4Patch(IndexType c) const = 0;
  virtual std::set<IndexType> CellNeighbours(IndexType i) const
  {
    std::cerr << "no CellNeighbours";
    abort();
    return std::set<IndexType>();
  }

  virtual int GetBoundaryCellOfCurved(IndexType iq) const { return -1; }

  virtual void AddShape(IndexType col, BoundaryFunction<2>* f)
  {
    std::cerr << "\"HierarchicalMesh::AddShape\" not written!" << std::endl;
    abort();
  }
  virtual void AddShape(IndexType col, BoundaryFunction<3>* f)
  {
    std::cerr << "\"HierarchicalMesh::AddShape\" not written!" << std::endl;
    abort();
  }

  void ShowOutput(IndexType i) const { _i_showoutput = i; }
  virtual void ProjectBoundary()
  {
    std::cout << "ProjectBoundary not written" << std::endl;
  }
};
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif
