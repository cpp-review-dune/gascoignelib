// WARNING DO NOT INCLUDE DIRECTLY

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

#include "../Common/vertex.h"
#include "../Interface/gascoigne.h"

#include "p4estdofhandlerbase.h"

namespace Gascoigne {

/// @brief
/// @tparam DEGREE of nodes
template<IndexType DEG>
class P4estDofHandler : public P4estDofHandlerBase
{
private:
  p4est_t* p4est; //< Not owned
  p4est_lnodes_t* lnodes;
  std::vector<std::array<MatrixEntryType, DIMENSION>>
    lnode_pos; //< Position of the lnodes

//   std::map<std::array<IndexType, DEG>, std::array<IndexType, DEG + 1>>
//     hn_edges;

// #if DIMENSION == 3
//   /// @brief Structur for hanging nodes in faces
//   std::map<std::array<IndexType, (DEG) * (DEG)>,
//            std::array<IndexType, (DEG + 1) * (DEG + 1)>>
//     hn_faces;
// #endif

  std::vector<std::array<MatrixEntryType, DIMENSION>> generate_lnode_pos(const ColumnStencil& _hn);
  void reorder_hanging_nodes();
  void insert_hn(std::vector<std::vector<IndexType>>& hn,
                 IndexType global_quad_id,
                 IndexType c1,
                 IndexType c2,
                 IndexType where);

public:
  P4estDofHandler(p4est_t* pforest);
  virtual ~P4estDofHandler();

  IndexType dimension() const;
  IndexType degree() const; //< degree of the lnodes

  IndexVector get_nodes_of_cell(IndexType cell) const;
  IndexType get_node_of_cell(IndexType cell, IndexType i) const;

  IndexType nnodes() const;
  IndexType nnothanging() const;
  IndexType nhanging() const;

  std::array<MatrixEntryType, DIMENSION> vertex(IndexType node_index) const;

  void write_vtk(std::string file_name,
                 double time,
                 GhostVectorAgent& gva,
                 std::vector<std::string> vectors) const;
};
}