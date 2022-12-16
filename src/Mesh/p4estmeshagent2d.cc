/**
 *
 * Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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

#ifndef P4_TO_P8
#include "p4estmeshagent2d.h"
#endif

#include <algorithm>
#include <iostream>

namespace Gascoigne {

void
p4est_quadrant_init_fn(p4est_t* p4est,
                       p4est_topidx_t which_tree,
                       p4est_quadrant_t* quadrant)
{}

/***
 * @param gridname File name of grid
 * @param prerefine Number of global refinments after initial creation
 * @param degree The Degree of the Finite elements
 */
P4estMeshAgent2d::P4estMeshAgent2d(const std::string& gridname,
                                   IndexType prerefine,
                                   IndexType degree)
{
  conn = p4est_connectivity_read_inp(gridname.c_str());
  if (conn == NULL) {
    P4EST_LERRORF("Failed to read a valid connectivity from %s\n",
                  gridname.c_str());
    sc_abort();
  }

  /* Create a forest that is not refined; it consists of the root octant. */
  pforest = p4est_new(sc_MPI_COMM_NULL,
                      conn,
                      sizeof(pquadrant_data_t),
                      p4est_quadrant_init_fn,
                      &pforest_data);

  global_refine(prerefine);

  /* Create the ghost layer to learn about parallel neighbors. */
  p4est_ghost_t* ghost = p4est_ghost_new(pforest, P4EST_CONNECT_FULL);

  /* Create a node numbering for continuous linear finite elements. */
  plnodes = p4est_lnodes_new(pforest, ghost, degree);

  /* Destroy the ghost structure -- no longer needed after node creation. */
  p4est_ghost_destroy(ghost);
  ghost = NULL;
}

P4estMeshAgent2d::~P4estMeshAgent2d()
{
  /* Destroy the p4est and the connectivity structure. */
  p4est_destroy(pforest);
  p4est_connectivity_destroy(conn);
  p4est_lnodes_destroy(plnodes);
}

/**
 * @return IndexType Number of quads/hexes in mesh
 */
IndexType
P4estMeshAgent2d::num_cells() const
{
  IndexType quad_count = 0;
  for (IndexType i = 0; i < pforest->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(pforest->trees, i);
    quad_count += tree->quadrants.elem_count;
  }
  return quad_count;
}

/**
 * @return the number of lnodes in the mesh
 */
IndexType
P4estMeshAgent2d::num_nodes() const
{
  return plnodes->num_local_elements;
}

/**
 * @param cell Index of the cell for witch the lnode index is looked up
 * @return IndexVector the lnode indices of the cell
 */
IndexVector
P4estMeshAgent2d::get_nodes_of_cell(IndexType cell) const
{
  IndexVector nodes(P4EST_CHILDREN);
  for (IndexType i = 0; i < P4EST_CHILDREN; ++i) {
    nodes[i] = plnodes->element_nodes[P4EST_CHILDREN * cell + i];
  }
  return nodes;
}

/**
 * @brief Generates a vtk output file
 *
 * @param fname Filename/path of vtk file to write
 */
void
P4estMeshAgent2d::write_vtk(const std::string& fname) const
{
  /* Write the forest to disk for visualization, one file per processor. */
  p4est_vtk_write_file(pforest, NULL, fname.c_str());
}

/**
 * @brief subdevides all cells of the mesh refine_level of times
 *
 * @param refine_level Number of subdevisions
 */
void
P4estMeshAgent2d::global_refine(IndexType refine_level)
{
  /** Callback function to decide on refinement.
   *
   * Refinement and coarsening is controlled by callback functions.
   * This function is called for every processor-local quadrant in order; its
   * return value is understood as a boolean refinement flag.
   *
   * Here we use uniform refinement.  Note that this function is not suitable
   * for recursive refinement and must be used in an iterative fashion.
   */
  auto refine_fn = [](p4est_t* p4est,
                      p4est_topidx_t which_tree,
                      p4est_quadrant_t* quadrant) -> int { return 1; };
  /* Refine the forest iteratively, load balancing at each iteration.
   * This is important when starting with an unrefined forest */

  for (IndexType level = 0; level < refine_level; ++level) {
    p4est_refine(pforest, 0, refine_fn, p4est_quadrant_init_fn);
  }
}

void
P4estMeshAgent2d::refine_cells(IndexVector& ref)
{
  for (IndexType i = 0; i < pforest->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(pforest->trees, i);
    for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
      p4est_quadrant_t* quadrant =
        p4est_quadrant_array_index(&(tree->quadrants), j);
      pquadrant_data_t* data = ((pquadrant_data_t*)quadrant->p.user_data);
      if (std::find(ref.begin(), ref.end(), tree->quadrants_offset + j) !=
          ref.end()) {
        data->refine = true;
      } else {
        data->refine = false;
      }
    }
  }

  auto refine_fn = [](p4est_t* p4est,
                      p4est_topidx_t which_tree,
                      p4est_quadrant_t* quadrant) -> int {
    if (((pquadrant_data_t*)quadrant->p.user_data)->refine == true) {
      ((pquadrant_data_t*)quadrant->p.user_data)->refine = false;
      return 1;
    }
    return 0;
  };

  p4est_refine(pforest, 0, refine_fn, p4est_quadrant_init_fn);
}

} // namespace Gascoigne