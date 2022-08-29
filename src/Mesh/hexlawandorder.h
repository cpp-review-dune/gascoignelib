/**
 *
 * Copyright (C) 2004 by the Gascoigne 3D authors
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

#ifndef __hexlawandorder_h
#define __hexlawandorder_h

#include "hex.h"
#include <map>

/*            2

  3__________6__________2
  |          |          |
  |          |          |
  |    3     |    2     |
  |          |          |
  |          |          |
3   7__________8__________5     1
  |          |          |
  |          |          |
  |    0     |    1     |
  |          |          |
  |          |          |
  0__________4__________1

             0
*/

/*---------------------------------------------------*/

namespace Gascoigne {
class HexLawAndOrder
{
protected:
  /* reference */
  std::vector<Hex>& hexs;

  /* data */
  std::array<EdgeVector, 12> lve, vice, childs_edge;
  std::array<FaceVector, 6> lvf, vicf, childs_face;
  std::array<IndexType, 8> cell_midpoint;
  std::array<std::array<IndexType, 3>, 8> ieoc;
  std::array<std::array<IndexType, 2>, 12> coif, lcfif;
  std::array<std::array<IndexType, 4>, 6> edgeofface;
  std::array<std::array<IndexType, 9>, 4> hnpf;

  void local_edge_index(EdgeVector& index, IndexType) const;
  void local_face_index(FaceVector& index, IndexType) const;
  void GetGlobalOuterFaceOfChild(FaceVector&,
                                 const Hex& f,
                                 IndexType c,
                                 IndexType e) const;

public:
  HexLawAndOrder(std::vector<Hex>&);

  void fill_corner_vertex_in_childs(const Hex&) const;
  void fill_middle_vertex_in_childs(const Hex&, IndexType) const;
  void fill_face_vertex_in_childs(const Hex&, IndexType, IndexType) const;
  void fill_edge_vertex_in_childs(const Hex&, IndexType, IndexType) const;

  IndexType face_vertex(const Hex&, IndexType) const;
  IndexType edge_vertex(const Hex&, IndexType) const;
  IndexType middle_vertex(const Hex&) const;
  void globalvertices_of_edge(const Hex& q,
                              std::array<IndexType, 2>& f,
                              IndexType ie) const;
  void globalvertices_of_face(const Hex& q,
                              std::array<IndexType, 4>& f,
                              IndexType ie) const;
  void LoadEdgeVerticesOfFace(const Hex& f,
                              IndexType face,
                              FaceVector& dst) const;
  void LoadFaceVertices(const Hex& f, std::array<IndexType, 6>& dst) const;

  // faces

  IndexType LocalChildFaceOfInnerFace(IndexType e, IndexType i) const
  {
    return lcfif[e][i];
  }
  IndexType ChildOfInnerFace(IndexType e, IndexType i) const
  {
    return coif[e][i];
  }
  IndexType ChildFace(IndexType e) const { return e; }
  IndexType ChildsOfFace(IndexType f, IndexType i) const
  {
    return childs_face[f][i];
  }
  IndexType local_face(const Hex&, const FaceVector&) const;
  IndexType GlobalInnerFace(IndexType c, IndexType i) const;
  IndexType GlobalChildFace(const FaceVector& edge,
                            IndexType q,
                            IndexType j) const;
  void GetFace(FaceVector&, IndexType h, IndexType e) const;
  void global_face_unsorted(FaceVector&, const Hex&, IndexType) const;
  void childs_of_face(FaceVector&, const Hex&, IndexType) const;
  void childs_of_global_face(FaceVector& child,
                             const Hex&,
                             const FaceVector&) const;

  std::pair<IndexType, IndexType> GetChildFaces(const FaceVector& bigedge,
                                                IndexType bighex,
                                                IndexType i) const;

  // edges

  void global_edge_unsorted(EdgeVector&, const Hex&, IndexType) const;
  IndexType GlobalChildEdge(const EdgeVector& edge,
                            IndexType q,
                            IndexType j) const;
  IndexType InnerEdge(const Hex&, IndexType i) const;
  IndexType InnerEdgeOfChild(IndexType c, IndexType i) const
  {
    return ieoc[c][i];
  }

  void load_face(FaceVector&, const Hex&, IndexType) const;
  IndexType local_face_index(IndexType, const FaceVector&) const;
  void globalfacechildren_of_father(std::vector<FaceVector>& faces,
                                    const Hex& f) const;
  IndexType LoadEdgeOfFace(const Hex& q,
                           const FaceVector& F,
                           IndexType e,
                           EdgeVector& E) const;
  IndexType LocalEdgeOfLocalFace(IndexType face, IndexType e) const
  {
    return edgeofface[face][e];
  }
  IndexType TestFaceOfOneChild(const Hex& f, const FaceVector& F) const;
  IndexType GetVertexOfEdge(IndexType iq,
                            const std::array<IndexType, 2>& edge) const;
  IndexType EdgeVertexOfFace(const Hex& q,
                             const FaceVector& F,
                             IndexType e) const;
  std::array<IndexType, 9> PatchVerticesOfFace(IndexType hex,
                                               IndexType face) const;
  std::array<IndexType, 9> GiveOrdering(const std::array<IndexType, 9>& F,
                                        const Hex& qfn) const;
};
} // namespace Gascoigne

#endif
