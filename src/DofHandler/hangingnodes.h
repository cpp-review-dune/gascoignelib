/*----------------------------   hangingnodes.h ---------------------------*/
/*      $Id:$                 */
#ifndef __hangingnodes_H
#define __hangingnodes_H
/*----------------------------   hangingnodes.h ---------------------------*/

/**
 *
 * Copyright (C) 2018 by the Gascoigne 3D authors
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

#include "../Common/entrymatrix.h"
#include "../Interface/hnstructureinterface.h"

/*-----------------------------------------*/

namespace Gascoigne {
// DIM:   spatial dimension
// NODES: nodes involved (Q1: 2, Q2: 3)

template<int DIM, int NODES>
class HangingNodes : public HNStructureInterface
{
  typedef std::array<IndexType, 3> EdgeVector;
  typedef std::array<IndexType, 9> FaceVector;

protected:
  const std::map<IndexType, EdgeVector>* edges;
  const std::map<IndexType, FaceVector>* faces;

  std::array<double, NODES> wei;
  std::array<double, NODES * NODES> fwei;

  std::array<std::array<IndexType, 3>, (DIM == 2) ? 4 : 12> lnoe;
  std::array<std::array<IndexType, (DIM == 2) ? 3 : 5>, (DIM == 2) ? 4 : 6>
    lnop;

public:
  HangingNodes();
  ~HangingNodes() {}

  void ReInit(const GascoigneMesh* M)
  {
    edges = M->GetHangingIndexHandler().GetStructure();
    assert(edges);
    if (DIM == 3) {
      faces = M->GetHangingIndexHandler().GetStructureFace();
      assert(faces);
    }
  }

  virtual void MatrixDiag(IndexType ncomp, MatrixInterface& A) const
  {
    nmatrix<double> M(ncomp);
    M.identity();

    for (auto it : *edges)
      A.entry_diag(it.first, M);
    if (DIM == 3)
      for (auto it : *faces)
        A.entry_diag(it.first, M);
  }

  virtual void SparseStructureDiag(SparseStructure* S) const
  {
    for (auto it : *edges)
      S->build_add(it.first, it.first);
    if (DIM == 3)
      for (auto it : *faces)
        S->build_add(it.first, it.first);
  }

  ///////////////////////////// Hanging nodes operations on vectors
  virtual void Average(GlobalVector& u) const
  {
    for (auto it : *edges) {
      u.zero_node(it.first);
      for (IndexType j = 0; j < NODES; ++j)
        u.add_node(it.first, wei[j], it.second[j]);
    }
    if (DIM == 3)
      for (auto it : *faces) {
        u.zero_node(it.first);
        for (IndexType j = 0; j < NODES * NODES; ++j)
          u.add_node(it.first, fwei[j], it.second[j]);
      }
  }

  void Distribute(GlobalVector& u) const
  {
    for (auto it : *edges) {
      for (IndexType j = 0; j < NODES; ++j)
        u.add_node(it.second[j], wei[j], it.first);
      u.zero_node(it.first);
    }
    if (DIM == 3)
      for (auto it : *faces) {
        for (IndexType j = 0; j < NODES * NODES; ++j)
          u.add_node(it.second[j], fwei[j], it.first);
        u.zero_node(it.first);
      }
  }

  void Zero(GlobalVector& u) const
  {
    for (auto it : *edges)
      u.zero_node(it.first);
    if (DIM == 3)
      for (auto it : *faces)
        u.zero_node(it.first);
  }

  bool ZeroCheck(const GlobalVector& u) const
  {
    for (auto it : *edges)
      for (ShortIndexType c = 0; c < u.ncomp(); c++)
        if (u(it.first, c) != 0.)
          return 1;
    if (DIM == 3)
      for (auto it : *faces)
        for (ShortIndexType c = 0; c < u.ncomp(); c++)
          if (u(it.first, c) != 0.)
            return 1;
    return 0;
  }

  IndexType nhnodes() const
  {
    if (DIM == 2)
      return edges->size();
    else
      return edges->size() + faces->size();
    ;
  }

  IndexType hanging(IndexType i) const
  {
    // hanging on line: returns 1 in 2d and 2 in 3d
    if (edges->find(i) != edges->end())
      return DIM - 1;
    if (DIM == 3)
      // hanging on face: returns 4
      if (faces->find(i) != faces->end())
        return 4;
    return 0;
  }

  std::array<IndexType, 2> GetHangingEdge(IndexType i) const
  {
    auto p = edges->find(i);
    assert(p != edges->end());

    std::array<IndexType, 2> Edge;
    for (IndexType j = 0; j < 2; ++j)
      Edge[j] = p->second[j];
    return Edge;
  }
  std::array<IndexType, 4> GetHangingFace(IndexType i) const
  {
    auto p = faces->find(i);
    assert(p != faces->end());

    std::array<IndexType, 4> Face;
    Face[0] = p->second[0];
    Face[1] = p->second[1];
    Face[2] = p->second[3];
    Face[3] = p->second[4];

    return Face;
  }
  const EdgeVector& regular_nodes(IndexType i) const
  {
    auto p = edges->find(i);
    assert(p != edges->end());
    return p->second;
  }

  void CondenseHanging(IndexVector& indices) const
  {
    if (DIM == 2)
      assert(indices.size() == NODES * NODES);
    else
      assert(indices.size() == NODES * NODES * NODES);
    // Q1
    if (NODES == 2) {
      // edges
      for (IndexType i = 0; i < indices.size(); ++i) {
        IndexType node = indices[i];
        auto edge = edges->find(node);
        if (edge != edges->end()) {
          // find node not in element
          for (IndexType j = 0; j < 2; ++j) {
            IndexType k = 0;
            for (k = 0; k < indices.size(); ++k)
              if (indices[k] == edge->second[j])
                break;
            if (k < indices.size()) // node found in element
              continue;             //
            indices[edge->first] = indices[edge->second[j]];
            break;
          }
        }
      }
      // faces
      if (DIM == 3)
        for (IndexType i = 0; i < indices.size(); ++i) {
          IndexType node = indices[i];
          auto face = faces->find(node);
          if (face != faces->end()) {
            // find node not in element
            for (IndexType j = 0; j < 4; ++j) {
              IndexType k = 0;
              for (k = 0; k < indices.size(); ++k)
                if (indices[k] == face->second[j])
                  break;
              if (k < indices.size()) // node found in element
                continue;             //
              indices[face->first] = indices[face->second[j]];
              break;
            }
          }
        }
    } else // NODES = 3
    {

      for (IndexType i = 0; i < indices.size(); ++i) {
        IndexType node = indices[i];
        auto edge = edges->find(node);
        if (edge != edges->end()) // node hangs, replace index i
          indices[i] = edge->second[NODES - 1];

        if (DIM == 3) {
          auto face = faces->find(node);
          if (face != faces->end()) // node hangs, replace index i
            indices[i] = face->second[NODES * NODES - 1];
        }
      }
    }
  }

  void CondenseHanging(EntryMatrix& E, IndexVector& indices) const;

  void CondenseHangingPatch(EntryMatrix& /*E*/, IndexVector& /*indices*/) const
  {
    //      assert(0);
  }
};
} // namespace Gascoigne

/*----------------------------   hangingnodes.h ---------------------------*/
/* end of #ifndef __hangingnodes_H */
#endif
/*----------------------------   hangingnodes.h ---------------------------*/
