/**
 *
 * Copyright (C) 2019,2020 by the Gascoigne 3D authors
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

#include "cgdofhandler.h"

#include "../Discretization/Q1/baseq12d.h"
#include "../Discretization/Q1/baseq13d.h"
#include "dofhandler.h"

/*-----------------------------------------*/

namespace Gascoigne {
/**
 * Structure to collect dofs on edges
 *
 * an edge is stored as two dof's (i1,i2) with i1<i2
 *
 **/
class PEdge
{

public:
  std::array<IndexType, 2> I;
  PEdge(IndexType a, IndexType b)
  {
    assert(a != b);
    I[0] = (a < b) ? a : b;
    I[1] = (a < b) ? b : a;
  }
  bool orient(IndexType a, IndexType b) const
  {
    if ((I[0] == a) && (I[1] == b))
      return true;
    else if ((I[0] == b) && (I[1] == a))
      return false;
    else {
      std::cerr << "edges do not match" << std::endl;
      abort();
    }
  }
  bool operator<(const PEdge x) const // lexicographic ordering of the edges
  {
    if (I[0] == x.I[0])
      return I[1] < x.I[1];
    return I[0] < x.I[0];
  }
};

class PrepareEdges
{
  const IndexType edge3d[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
                                    { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
                                    { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };

public:
  std::map<PEdge, IndexType> edges;

  IndexType nedges() const { return edges.size(); }

  template<int DIM>
  void InitEdges(const DofHandler<DIM>& DH)
  {
    const IndexVector& nc = DH.GetCellVector();
    IndexType nedges = 0;
    for (IndexType c = 0; c < DH.ncells(); ++c) {
      IndexType nd = (DIM == 2) ? 4 : 8;
      IndexType start = c * nd;
      if (DIM == 2)
        for (IndexType e = 0; e < 4; ++e) {
          PEdge edge(nc[start + e], nc[start + (e + 1) % nd]);
          if (edges.find(edge) == edges.end()) {
            edges[edge] = nedges;
            ++nedges;
          }
        }
      else if (DIM == 3)
        for (IndexType e = 0; e < 12; ++e) {
          PEdge edge(nc[start + edge3d[e][0]], nc[start + edge3d[e][1]]);
          if (edges.find(edge) == edges.end()) {
            edges[edge] = nedges;
            ++nedges;
          }
        }
    }
  }

  IndexType dof_on_edge(IndexType i1,
                        IndexType i2,
                        IndexType j,
                        IndexType M) const
  {
    assert(j >= 0);
    assert(j < M - 2);
    PEdge edge(i1, i2);
    auto it = edges.find(edge);
    assert(it != edges.end());
    IndexType edgeid = it->second;
    bool orientation = edge.orient(i1, i2);
    if (orientation)
      return edgeid * (M - 2) + j;
    else
      return edgeid * (M - 2) + (M - 3) - j;
  }
};

/**
 * Structure to collect dofs on faces
 *
 * A face is stored as
 *
 * i3  i4
 *
 * i1  i2
 *
 * with i1 being the smallest dof index and i2<i3
 *
 **/

class PFace
{
public:
  std::array<IndexType, 4> I;
  PFace(IndexType a, IndexType b, IndexType c, IndexType d)
  {
    if ((a < b) && (a < c) && (a < d))
      I = { a, b, c, d };
    else if ((b < a) && (b < c) && (b < d))
      I = { b, d, a, c };
    else if ((c < a) && (c < b) && (c < d))
      I = { c, a, d, b };
    else if ((d < a) && (d < b) && (d < c))
      I = { d, c, b, a };
    else
      abort();
    if (I[2] < I[1])
      std::swap(I[1], I[2]);
  }
  /**
   *  std numbering of face with. i0< i1,i2,i3 and i1<i2
   *  2--3
   *  |  |
   *  0--1
   *
   *  orientation
   *
   *   +1      +2     +3     +4
   *  +y      -x     -y     +x
   *  2--3    0--2   1--0   3--1     y
   *  |  |    |  |   |  |   |  |     |
   *  0--1+x  1--3+y 3--2-x 2--0-y   +--x
   *
   *   -1    -2     -3     -4
   *  +x      +y     -x     -y
   *  1--3    3--2   2--0   0--1     y
   *  |  |    |  |   |  |   |  |     |
   *  0--2+y  1--0-x 3--1-y 2--3+x   +--x
   *
   **/
  int orientation(IndexType a, IndexType b) const
  {
    if ((I[0] == a) && (I[1] == b))
      return +1;
    else if ((I[0] == a) && (I[2] == b))
      return -1;
    else if ((I[1] == a) && (I[3] == b))
      return +2;
    else if ((I[1] == a) && (I[0] == b))
      return -2;
    else if ((I[3] == a) && (I[2] == b))
      return +3;
    else if ((I[3] == a) && (I[1] == b))
      return -3;
    else if ((I[2] == a) && (I[0] == b))
      return +4;
    else if ((I[2] == a) && (I[3] == b))
      return -4;
    else
      abort();
  }
  bool operator<(const PFace& x) const // lexicographic ordering of the edges
  {
    if (I[0] == x.I[0]) {
      if (I[1] == x.I[1]) {
        if (I[2] == x.I[2])
          return I[3] < x.I[3];
        else
          return I[2] < x.I[2];
      } else
        return I[1] < x.I[1];
    } else
      return I[0] < x.I[0];
  }
};

class PrepareFaces
{
  /// four nodes each define one face.
  const IndexType face3d[6][4] = { { 0, 1, 3, 2 }, { 1, 5, 2, 6 },
                                   { 4, 5, 7, 6 }, { 0, 4, 3, 7 },
                                   { 0, 1, 4, 5 }, { 3, 2, 7, 6 } };

public:
  std::map<PFace, IndexType> faces;

  IndexType nfaces() const { return faces.size(); }

  template<int DIM>
  void InitFaces(const DofHandler<DIM>& DH)
  {
    assert(DIM == 3);
    const IndexVector& nc = DH.GetCellVector();
    assert(nc.size() == 8 * DH.ncells());

    // DIM==3
    IndexType nfaces = 0;
    for (IndexType c = 0; c < DH.ncells(); ++c) {
      IndexType nd = 8;         // dofs per element
      IndexType start = c * nd; // first index of dof
      for (IndexType f = 0; f < 6; ++f) {
        PFace face(nc[start + face3d[f][0]],
                   nc[start + face3d[f][1]],
                   nc[start + face3d[f][2]],
                   nc[start + face3d[f][3]]);
        if (faces.find(face) == faces.end()) {
          faces[face] = nfaces;
          ++nfaces;
        }
      }
    }
  }
  /**
   *  std numbering of face with. i0< i1,i2,i3 and i1<i2
   *  2--3
   *  |  |
   *  0--1
   *
   *  orientation
   *
   *   +1      +2     +3     +4
   *  +y      -x     -y     +x
   *  2--3    0--2   1--0   3--1     y
   *  |  |    |  |   |  |   |  |     |
   *  0--1+x  1--3+y 3--2-x 2--0-y   +--x
   *
   *   -1    -2     -3     -4
   *  +x      +y     -x     -y
   *  1--3    3--2   2--0   0--1     y
   *  |  |    |  |   |  |   |  |     |
   *  0--2+y  1--0-x 3--1-y 2--3+x   +--x
   *
   **/
  IndexType dof_on_face(IndexType a,
                        IndexType b,
                        IndexType c,
                        IndexType d,
                        IndexType ix,
                        IndexType iy,
                        IndexType M) const
  {
    assert(ix >= 0);
    assert(iy >= 0);
    assert(ix <= M - 2);
    assert(iy <= M - 2);
    PFace face(a, b, c, d);
    auto it = faces.find(face);
    assert(it != faces.end());
    IndexType faceid = it->second;
    int o = face.orientation(a, b);
    if (o == 1)
      return faceid * (M - 2) * (M - 2) + (M - 2) * iy + ix;
    else if (o == 2)
      return faceid * (M - 2) * (M - 2) + (M - 2) * (M - 3 - ix) + iy;
    else if (o == 3)
      return faceid * (M - 2) * (M - 2) + (M - 2) * (M - 3 - iy) + (M - 3 - ix);
    else if (o == 4)
      return faceid * (M - 2) * (M - 2) + (M - 2) * ix + (M - 3 - iy);
    else if (o == -1)
      return faceid * (M - 2) * (M - 2) + (M - 2) * ix + iy;
    else if (o == -2)
      return faceid * (M - 2) * (M - 2) + (M - 2) * iy + (M - 3 - ix);
    else if (o == -3)
      return faceid * (M - 2) * (M - 2) + (M - 2) * (M - 3 - ix) + (M - 3 - iy);
    else if (o == -4)
      return faceid * (M - 2) * (M - 2) + (M - 2) * (M - 3 - iy) + ix;
    else
      abort();
  }
};

/**
 *
 *
 *
 *
 **/
template<int DIM, int M>
void
CGDofHandler<DIM, M>::InitFromGascoigneMesh(const DofHandler<DIM>& GM)
{
  assert((DIM == 2) || (DIM == 3));
  assert(M >= 2);

  const IndexVector& GM_nc = GM.GetCellVector();
  const std::vector<Vertex<DIM>>& GM_nv = GM.GetVertexVector();
  IndexType GM_ncells = GM.ncells();
  IndexType GM_ndofs = GM_nv.size();

  ////////// preparation
  // first, dof's on the nodes. This is the Gascoignemesh structure, e.g. first
  // GM_ndofs nodes then, dof's on the edges - only for M>2 then, dof's on faces
  // in 3D  - only for M>2 finally, dof's within the elements - only for M>2

  // edge list
  PrepareEdges prepedges;
  prepedges.InitEdges(GM);

  PrepareFaces prepfaces;
  if (DIM == 3)
    prepfaces.InitFaces(GM);

  IndexType dofs_on_nodes = static_cast<IndexType>(GM_ndofs);
  IndexType dofs_on_edges =
    (M - 2) * static_cast<IndexType>(prepedges.nedges());
  IndexType dofs_on_faces = 0;
  if (DIM == 3)
    dofs_on_faces =
      (M - 2) * (M - 2) * static_cast<IndexType>(prepfaces.nfaces());

  IndexType dofs_on_elements =
    static_cast<IndexType>(pow(M - 2, DIM)) * static_cast<IndexType>(GM_ncells);

  IndexType ndofs =
    dofs_on_nodes + dofs_on_edges + dofs_on_faces + dofs_on_elements;
  std::cout << "CGDofHandler " << ndofs << "\t" << dofs_on_nodes << " "
            << dofs_on_edges << " " << dofs_on_faces << " " << dofs_on_elements
            << std::endl;

  // create list of DOF's
  nc.resize(dofs_per_element() * GM_ncells, -1);

  // CGDofHandler - Element == GascoigneMesh - Cell
  for (IndexType element = 0; element < GM_ncells; ++element) {
    IndexType start = element * dofs_per_element();
    IndexType GM_start = element * ((DIM == 2) ? 4 : 8);

    if (DIM == 2)
      for (IndexType iy = 0; iy < M; ++iy)
        for (IndexType ix = 0; ix < M; ++ix) {
          IndexType dof = start + iy * M + ix;
          // dofs in the 4 nodes
          if ((ix == 0) && (iy == 0))
            nc[dof] = GM_nc[GM_start + 0]; // ll
          else if ((ix == (M - 1)) && (iy == 0))
            nc[dof] = GM_nc[GM_start + 1]; // lr
          else if ((ix == 0) && (iy == (M - 1)))
            nc[dof] = GM_nc[GM_start + 3]; // ul
          else if ((ix == (M - 1)) && (iy == (M - 1)))
            nc[dof] = GM_nc[GM_start + 2]; // ur
          // dofs on the four lines
          else if (iy == 0) // lower
            nc[dof] = dofs_on_nodes +
                      prepedges.dof_on_edge(
                        GM_nc[GM_start], GM_nc[GM_start + 1], ix - 1, M);
          else if (ix == (M - 1)) // right
            nc[dof] = dofs_on_nodes +
                      prepedges.dof_on_edge(
                        GM_nc[GM_start + 1], GM_nc[GM_start + 2], iy - 1, M);
          else if (iy == (M - 1)) // top
            nc[dof] = dofs_on_nodes +
                      prepedges.dof_on_edge(
                        GM_nc[GM_start + 3], GM_nc[GM_start + 2], ix - 1, M);
          else if (ix == 0) // left
            nc[dof] = dofs_on_nodes +
                      prepedges.dof_on_edge(
                        GM_nc[GM_start], GM_nc[GM_start + 3], iy - 1, M);
          // nodes within the domain
          else
            nc[dof] = dofs_on_nodes + dofs_on_edges +
                      (M - 2) * (M - 2) * element + (M - 2) * (iy - 1) +
                      (ix - 1);
        }
    else if (DIM == 3)
      /**
       *
       *    7  ---   6
       *   /|       /|
       *  3  ----  2 |    y   z
       *  | |      | |    |  /
       *  | 4  ----| 5    | /
       *  |/       |/     |/
       *  0  ----  1      +---- x
       *
       * 8 node - dofs
       * 12 edges (0,1), (4,5), (3,2), (7,6)
       *          (0,4), (3,7), (1,5), (2,6)
       *          (0,3), (1,2), (4,7), (5,6)
       * 6 faces  (0,1,3,2), (4,5,7,6)
       *          (0,4,3,7), (1,5,2,6)
       *          (0,1,4,5), (3,2,7,6)
       **/
      for (IndexType iz = 0; iz < M; ++iz)
        for (IndexType iy = 0; iy < M; ++iy)
          for (IndexType ix = 0; ix < M; ++ix) {
            IndexType dof = start + iz * M * M + iy * M + ix; // index of dof

            // dofs in the 8 nodes
            if ((ix == 0) && (iy == 0) && (iz == 0))
              nc[dof] = GM_nc[GM_start + 0]; // lll
            else if ((ix == M - 1) && (iy == 0) && (iz == 0))
              nc[dof] = GM_nc[GM_start + 1]; // rll
            else if ((ix == 0) && (iy == M - 1) && (iz == 0))
              nc[dof] = GM_nc[GM_start + 3]; // lrl
            else if ((ix == M - 1) && (iy == M - 1) && (iz == 0))
              nc[dof] = GM_nc[GM_start + 2]; // rrl
            else if ((ix == 0) && (iy == 0) && (iz == M - 1))
              nc[dof] = GM_nc[GM_start + 4]; // llr
            else if ((ix == M - 1) && (iy == 0) && (iz == M - 1))
              nc[dof] = GM_nc[GM_start + 5]; // rlr
            else if ((ix == 0) && (iy == M - 1) && (iz == M - 1))
              nc[dof] = GM_nc[GM_start + 7]; // lrr
            else if ((ix == M - 1) && (iy == M - 1) && (iz == M - 1))
              nc[dof] = GM_nc[GM_start + 6]; // rrr
            // edges
            else if ((iy == 0) && (iz == 0))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 0], GM_nc[GM_start + 1], ix - 1, M);
            else if ((iy == 0) && (iz == M - 1))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 4], GM_nc[GM_start + 5], ix - 1, M);
            else if ((iy == M - 1) && (iz == 0))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 3], GM_nc[GM_start + 2], ix - 1, M);
            else if ((iy == M - 1) && (iz == M - 1))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 7], GM_nc[GM_start + 6], ix - 1, M);
            else if ((ix == 0) && (iy == 0))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 0], GM_nc[GM_start + 4], iz - 1, M);
            else if ((ix == 0) && (iy == M - 1))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 3], GM_nc[GM_start + 7], iz - 1, M);
            else if ((ix == M - 1) && (iy == 0))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 1], GM_nc[GM_start + 5], iz - 1, M);
            else if ((ix == M - 1) && (iy == M - 1))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 2], GM_nc[GM_start + 6], iz - 1, M);
            else if ((iz == 0) && (ix == 0))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 0], GM_nc[GM_start + 3], iy - 1, M);
            else if ((iz == 0) && (ix == M - 1))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 1], GM_nc[GM_start + 2], iy - 1, M);
            else if ((iz == M - 1) && (ix == 0))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 4], GM_nc[GM_start + 7], iy - 1, M);
            else if ((iz == M - 1) && (ix == M - 1))
              nc[dof] = dofs_on_nodes +
                        prepedges.dof_on_edge(
                          GM_nc[GM_start + 5], GM_nc[GM_start + 6], iy - 1, M);
            // faces
            else if (iz == 0)
              nc[dof] = dofs_on_nodes + dofs_on_edges +
                        prepfaces.dof_on_face(GM_nc[GM_start + 0],
                                              GM_nc[GM_start + 1],
                                              GM_nc[GM_start + 3],
                                              GM_nc[GM_start + 2],
                                              ix - 1,
                                              iy - 1,
                                              M);
            else if (iz == M - 1)
              nc[dof] = dofs_on_nodes + dofs_on_edges +
                        prepfaces.dof_on_face(GM_nc[GM_start + 4],
                                              GM_nc[GM_start + 5],
                                              GM_nc[GM_start + 7],
                                              GM_nc[GM_start + 6],
                                              ix - 1,
                                              iy - 1,
                                              M);
            else if (ix == 0)
              nc[dof] = dofs_on_nodes + dofs_on_edges +
                        prepfaces.dof_on_face(GM_nc[GM_start + 0],
                                              GM_nc[GM_start + 4],
                                              GM_nc[GM_start + 3],
                                              GM_nc[GM_start + 7],
                                              iz - 1,
                                              iy - 1,
                                              M);
            else if (ix == M - 1)
              nc[dof] = dofs_on_nodes + dofs_on_edges +
                        prepfaces.dof_on_face(GM_nc[GM_start + 1],
                                              GM_nc[GM_start + 5],
                                              GM_nc[GM_start + 2],
                                              GM_nc[GM_start + 6],
                                              iz - 1,
                                              iy - 1,
                                              M);
            else if (iy == 0)
              nc[dof] = dofs_on_nodes + dofs_on_edges +
                        prepfaces.dof_on_face(GM_nc[GM_start + 0],
                                              GM_nc[GM_start + 1],
                                              GM_nc[GM_start + 4],
                                              GM_nc[GM_start + 5],
                                              ix - 1,
                                              iz - 1,
                                              M);
            else if (iy == M - 1)
              nc[dof] = dofs_on_nodes + dofs_on_edges +
                        prepfaces.dof_on_face(GM_nc[GM_start + 3],
                                              GM_nc[GM_start + 2],
                                              GM_nc[GM_start + 7],
                                              GM_nc[GM_start + 6],
                                              ix - 1,
                                              iz - 1,
                                              M);
            else
              nc[dof] = dofs_on_nodes + dofs_on_edges + dofs_on_faces +
                        (M - 2) * (M - 2) * (M - 2) * element +
                        (M - 2) * (M - 2) * (iz - 1) + (M - 2) * (iy - 1) +
                        (ix - 1);
          }
    else
      abort();
  }
  for (auto it : nc)
    assert(it >= 0);

  //////////////////// check that all dofs are found
  std::vector<bool> df(ndofs, false);
  for (auto it : nc)
    df[it] = true;
  for (auto it : df)
    assert(it);

  //////////////////// create vertex coordinates
  nx.resize(ndofs);
  if (DIM == 2) {
    BaseQ12d B;
    for (IndexType el = 0; el < GM_ncells; ++el) {
      IndexType GM_start = el * ((DIM == 2) ? 4 : 8);
      for (IndexType iy = 0; iy < M; ++iy)
        for (IndexType ix = 0; ix < M; ++ix) {
          Vertex<2> s;
          s.x() = static_cast<double>(ix) / static_cast<double>(M - 1);
          s.y() = static_cast<double>(iy) / static_cast<double>(M - 1);
          B.point(s);
          nx[nc[dofs_per_element() * el + iy * M + ix]].equ(
            B.phi(0),
            GM_nv[GM_nc[GM_start + 0]],
            B.phi(1),
            GM_nv[GM_nc[GM_start + 1]],
            B.phi(2),
            GM_nv[GM_nc[GM_start + 3]],
            B.phi(3),
            GM_nv[GM_nc[GM_start + 2]]);
        }
    }
  } else if (DIM == 3) {
    BaseQ13d B;
    for (IndexType el = 0; el < GM_ncells; ++el) {
      IndexType GM_start = el * ((DIM == 2) ? 4 : 8);
      for (IndexType iz = 0; iz < M; ++iz)
        for (IndexType iy = 0; iy < M; ++iy)
          for (IndexType ix = 0; ix < M; ++ix) {
            Vertex<3> s;
            s.x() = static_cast<double>(ix) / static_cast<double>(M - 1);
            s.y() = static_cast<double>(iy) / static_cast<double>(M - 1);
            s.z() = static_cast<double>(iz) / static_cast<double>(M - 1);
            B.point(s);
            nx[nc[dofs_per_element() * el + iz * M * M + iy * M + ix]].equ(
              { B.phi(0),
                B.phi(1),
                B.phi(2),
                B.phi(3),
                B.phi(4),
                B.phi(5),
                B.phi(6),
                B.phi(7) },
              { GM_nv[GM_nc[GM_start + 0]],
                GM_nv[GM_nc[GM_start + 1]],    //     7  ---   6
                GM_nv[GM_nc[GM_start + 3]],    //    /|       /|
                GM_nv[GM_nc[GM_start + 2]],    //   3  ----  2 |
                GM_nv[GM_nc[GM_start + 4]],    //   | |      | |
                GM_nv[GM_nc[GM_start + 5]],    //   | 4  ----| 5
                GM_nv[GM_nc[GM_start + 7]],    //   |/       |/
                GM_nv[GM_nc[GM_start + 6]] }); //   0  ----  1
          }
    }
  }

  // // DIAG: output all coordinates
  // for (IndexType el=0;el<GM_ncells;++el)
  //   {
  // 	if (DIM==2)
  // 	  for (IndexType iy=0;iy<M;++iy)
  // 	    for (IndexType ix=0;ix<M;++ix)
  // 	      std::cerr << nx[nc[dofs_per_element()*el+iy*M+ix]].x() << " "
  // 			<< nx[nc[dofs_per_element()*el+iy*M+ix]].y() <<
  // std::endl; 	else 	  for (IndexType iz=0;iz<M;++iz) 	    for
  // (IndexType iy=0;iy<M;++iy) 	      for (IndexType ix=0;ix<M;++ix)
  // std::cerr << nx[nc[dofs_per_element()*el+iz*M*M+iy*M+ix]].x() << " "
  // 			  << nx[nc[dofs_per_element()*el+iz*M*M+iy*M+ix]].y() <<
  // "
  // "
  // 			  << nx[nc[dofs_per_element()*el+iz*M*M+iy*M+ix]].z() <<
  // "
  // "
  // << std::endl; 	std::cerr << std::endl << std::endl;
  //   }

  ////////////////////////////////////////////////// BoundaryHandler
  BoundaryIndexHandler& BIH = GetBoundaryIndexHandler();
  const BoundaryIndexHandler& GM_BIH = GM.GetBoundaryIndexHandler();
  typedef std::map<IndexType, IndexVector> VecMap;

  IndexSet& colors = BIH.GetColors();
  VecMap& vertex = BIH.GetVertex();
  VecMap& cell = BIH.GetCell();
  VecMap& local = BIH.GetLocal();
  //    VecMap& patch = BIH.GetPatch();
  //    VecMap& localpatch = BIH.GetLocalPatch();

  const IndexSet& GM_colors = GM_BIH.GetColors();
  //    const VecMap& GM_vertex = GM_BIH.GetVertex();
  const VecMap& GM_cell = GM_BIH.GetCell();
  const VecMap& GM_local = GM_BIH.GetLocal();
  //    const VecMap& GM_patch  = GM_BIH.GetPatch();
  //    const VecMap& GM_localpatch = GM_BIH.GetLocalPatch();

  // Colors. Simple copy
  colors = GM_colors;
  // Cell & Local - Simple copy
  cell = GM_cell; // klappt das?
  assert(cell.size() == GM_cell.size());
  if (cell.size() > 0)
    assert(cell.begin()->second.size() == GM_cell.begin()->second.size());
  local = GM_local;
  // Vertices. New construction
  constexpr IndexType lvf[6][4] = {
    { 0, 1, 2, 3 }, { 1, 5, 6, 2 }, // 3d: nodes on face
    { 2, 6, 7, 3 }, { 3, 7, 4, 0 }, { 0, 4, 5, 1 }, { 4, 7, 6, 5 }
  };
  constexpr IndexType lve[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
                                     { 4, 5 }, { 5, 6 }, // 3d: nodes on edges
                                     { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 },
                                     { 2, 6 }, { 3, 7 } };

  constexpr IndexType eof[6][4] = { { 0, 1, 2, 3 },   { 9, 5, 10, 1 },
                                    { 10, 6, 11, 2 }, // 3d: edge of face
                                    { 11, 7, 8, 3 },  { 8, 4, 9, 0 },
                                    { 7, 6, 5, 4 } };

  for (auto col : colors) {
    const IndexVector& GM_cells = GM_BIH.Cells(col);
    const IndexVector& GM_local = GM_BIH.Localind(col);
    assert(GM_cells.size() == GM_local.size());
    std::set<IndexType> dofoncolor;
    for (IndexType i = 0; i < GM_cells.size(); ++i) {
      IndexType GM_c = GM_cells[i];
      IndexType GM_l = GM_local[i];
      IndexType GM_s =
        GM_c * ((DIM == 2) ? 4 : 8); // start index for nodes on edges

      if (DIM == 2) {
        //  (number of of GascoigneMesh)
        //  3 -- 2 -- 2
        //  |         |
        //  3         1
        //  |         |
        //  0 -- 0 -- 1

        // dof's in nodes
        dofoncolor.insert(GM_nc[GM_s + GM_l]);
        dofoncolor.insert(GM_nc[GM_s + (GM_l + 1) % 4]);
        // dofs along line
        for (IndexType m = 1; m < M - 1; ++m)
          dofoncolor.insert(
            dofs_on_nodes +
            prepedges.dof_on_edge(
              GM_nc[GM_s + GM_l], GM_nc[GM_s + (GM_l + 1) % 4], m - 1, M));
      } else {
        // nodes
        for (IndexType i = 0; i < 4; ++i)
          dofoncolor.insert(GM_nc[GM_s + lvf[GM_l][i]]);
        // edge
        for (IndexType i = 0; i < 4; ++i) {
          IndexType ei = eof[GM_l][i];
          PEdge edge(GM_nc[GM_s + lve[ei][0]], GM_nc[GM_s + lve[ei][1]]);
          auto ite = prepedges.edges.find(edge);
          assert(ite != prepedges.edges.end());
          for (IndexType ix = 0; ix < M - 2; ++ix)
            dofoncolor.insert(
              dofs_on_nodes +
              prepedges.dof_on_edge(
                GM_nc[GM_s + lve[ei][0]], GM_nc[GM_s + lve[ei][1]], ix, M));
        }
        // face
        PFace face(GM_nc[GM_s + lvf[GM_l][0]],
                   GM_nc[GM_s + lvf[GM_l][1]],
                   GM_nc[GM_s + lvf[GM_l][3]],
                   GM_nc[GM_s + lvf[GM_l][2]]);
        auto it = prepfaces.faces.find(face);
        assert(it != prepfaces.faces.end());
        for (IndexType iy = 0; iy < M - 2; ++iy)
          for (IndexType ix = 0; ix < M - 2; ++ix)
            dofoncolor.insert(dofs_on_nodes + dofs_on_edges +
                              prepfaces.dof_on_face(GM_nc[GM_s + lvf[GM_l][0]],
                                                    GM_nc[GM_s + lvf[GM_l][1]],
                                                    GM_nc[GM_s + lvf[GM_l][3]],
                                                    GM_nc[GM_s + lvf[GM_l][2]],
                                                    ix,
                                                    iy,
                                                    M));
      }
    }
    for (auto v : dofoncolor)
      vertex[col].push_back(v);
  }
  ////////////////////////////////////////////////// Material
  mat = GM.GetMaterialVector();
}
} // namespace Gascoigne

template class Gascoigne::CGDofHandler<2, 2>;
template class Gascoigne::CGDofHandler<2, 3>;
template class Gascoigne::CGDofHandler<2, 5>;

template class Gascoigne::CGDofHandler<3, 2>;
template class Gascoigne::CGDofHandler<3, 3>;
template class Gascoigne::CGDofHandler<3, 5>;
