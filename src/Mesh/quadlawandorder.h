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

#ifndef __quadlawandorder_h
#define __quadlawandorder_h

#include "cell.h"
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
class QuadLawAndOrder
{
protected:
  typedef std::map<IndexType, IndexType> LocVertexLocEdge;
  typedef std::vector<LocVertexLocEdge> LocVertexLocVertexLocEdge;
  typedef std::array<IndexType, 2> QuadVector;

  // Daten fuer Suchen von kindern an hang

  LocVertexLocVertexLocEdge lvlvle;

  // Referenzen fuer globale fkts

  std::vector<Quad>& quads;

  std::array<EdgeVector, 4> childs_edge, vice;
  std::array<IndexType, 4> child_point_cell, child_point_vertex;
  std::array<IndexType, 9> gc, gv;

  std::array<std::array<IndexType, 2>, 4> ieoc, oeoc;

  IndexType local_edge(const Quad& f, const EdgeVector& globaledge) const;

public:
  QuadLawAndOrder(std::vector<Quad>& q);

  IndexType cell_midpoint(IndexType i) const { return (i + 2) % 4; }
  IndexType global_index(const Quad& q, IndexType i) const;

  void fill_corner_vertex_in_childs(const Quad& f) const;
  void fill_edge_vertex_in_childs(const Quad& f,
                                  IndexType e,
                                  IndexType i) const;
  void fill_middle_vertex_in_childs(const Quad& f, IndexType i) const;

  //   edges

  IndexType ChildEdge(IndexType e) const { return e; }
  IndexType ChildsOfEdge(IndexType e, IndexType i) const
  {
    return childs_edge[e][i];
  }
  IndexType InnerEdgeOfChild(IndexType c, IndexType i) const
  {
    return ieoc[c][i];
  }
  IndexType OuterEdgeOfChild(IndexType c, IndexType i) const
  {
    return oeoc[c][i];
  }
  IndexType GlobalInnerEdge(IndexType c, IndexType i) const;

  std::pair<IndexType, IndexType> GetChildEdges(EdgeVector& edge,
                                                const EdgeVector& bigedge,
                                                IndexType hanging,
                                                IndexType bigquad,
                                                IndexType i) const;

  IndexType GlobalChildEdge(const EdgeVector& edge,
                            IndexType q,
                            IndexType j) const;
  void local_edge_index(EdgeVector& index, IndexType edge) const;
  IndexType local_edge_index(IndexType, const EdgeVector&) const;

  IndexType middle_vertex(const Quad& f) const;
  IndexType edge_vertex(const Quad& f, IndexType edge) const;

  /* for boundaries */
  void childs_of_edge(QuadVector& child, const Quad& f, IndexType edge) const;

  /* for regular */
  void childs_of_global_edge(QuadVector& child,
                             const Quad& f,
                             const EdgeVector& globaledge) const;

  void globaledgechildren_of_father(std::vector<EdgeVector>& edges,
                                    const Quad& f) const;

  /* for hierarchicalmesh / find in liehanglist */
  void global_edge_unsorted(std::array<IndexType, 2>& lineglob,
                            const Quad& q,
                            IndexType edge) const;

  /* fuer mginterpolator */
  void globalvertices_of_edge(const Quad&, EdgeVector&, IndexType) const;
};
} // namespace Gascoigne

#endif
