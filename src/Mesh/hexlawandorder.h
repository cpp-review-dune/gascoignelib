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
  /* typedef's */
  typedef std::array<int, 2> EdgeVector;
  typedef std::array<int, 4> FaceVector;

  /* reference */
  std::vector<Hex>& hexs;

  /* data */
  std::array<EdgeVector, 12> lve, vice, childs_edge;
  std::array<FaceVector, 6> lvf, vicf, childs_face;
  std::array<int, 8> cell_midpoint;
  std::array<std::array<int, 3>, 8> ieoc;
  std::array<std::array<int, 2>, 12> coif, lcfif;
  std::array<std::array<int, 4>, 6> edgeofface;
  std::array<std::array<int, 9>, 4> hnpf;

  void local_edge_index(EdgeVector& index, int) const;
  void local_face_index(FaceVector& index, int) const;
  void GetGlobalOuterFaceOfChild(FaceVector&, const Hex& f, int c, int e) const;

public:
  HexLawAndOrder(std::vector<Hex>&);

  void fill_corner_vertex_in_childs(const Hex&) const;
  void fill_middle_vertex_in_childs(const Hex&, int) const;
  void fill_face_vertex_in_childs(const Hex&, int, int) const;
  void fill_edge_vertex_in_childs(const Hex&, int, int) const;

  int face_vertex(const Hex&, int) const;
  int edge_vertex(const Hex&, int) const;
  int middle_vertex(const Hex&) const;
  void globalvertices_of_edge(const Hex& q,
                              std::array<int, 2>& f,
                              int ie) const;
  void globalvertices_of_face(const Hex& q,
                              std::array<int, 4>& f,
                              int ie) const;
  void LoadEdgeVerticesOfFace(const Hex& f, int face, FaceVector& dst) const;
  void LoadFaceVertices(const Hex& f, std::array<int, 6>& dst) const;

  // faces

  int LocalChildFaceOfInnerFace(int e, int i) const { return lcfif[e][i]; }
  int ChildOfInnerFace(int e, int i) const { return coif[e][i]; }
  int ChildFace(int e) const { return e; }
  int ChildsOfFace(int f, int i) const { return childs_face[f][i]; }
  int local_face(const Hex&, const FaceVector&) const;
  int GlobalInnerFace(int c, int i) const;
  int GlobalChildFace(const FaceVector& edge, int q, int j) const;
  void GetFace(FaceVector&, int h, int e) const;
  void global_face_unsorted(FaceVector&, const Hex&, int) const;
  void childs_of_face(FaceVector&, const Hex&, int) const;
  void childs_of_global_face(FaceVector& child,
                             const Hex&,
                             const FaceVector&) const;

  std::pair<int, int> GetChildFaces(const FaceVector& bigedge,
                                    int bighex,
                                    int i) const;

  // edges

  void global_edge_unsorted(EdgeVector&, const Hex&, int) const;
  int GlobalChildEdge(const EdgeVector& edge, int q, int j) const;
  int InnerEdge(const Hex&, int i) const;
  int InnerEdgeOfChild(int c, int i) const { return ieoc[c][i]; }

  void load_face(FaceVector&, const Hex&, int) const;
  int local_face_index(int, const FaceVector&) const;
  void globalfacechildren_of_father(std::vector<FaceVector>& faces,
                                    const Hex& f) const;
  int LoadEdgeOfFace(const Hex& q,
                     const FaceVector& F,
                     int e,
                     EdgeVector& E) const;
  int LocalEdgeOfLocalFace(int face, int e) const
  {
    return edgeofface[face][e];
  }
  int TestFaceOfOneChild(const Hex& f, const FaceVector& F) const;
  int GetVertexOfEdge(int iq, const std::array<int, 2>& edge) const;
  int EdgeVertexOfFace(const Hex& q, const FaceVector& F, int e) const;
  std::array<int, 9> PatchVerticesOfFace(int hex, int face) const;
  std::array<int, 9> GiveOrdering(const std::array<int, 9>& F,
                                  const Hex& qfn) const;
};
} // namespace Gascoigne

#endif
