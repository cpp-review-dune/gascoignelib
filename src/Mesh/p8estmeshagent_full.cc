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

#include "pforestmeshagent.h"

#include <algorithm>
#include <iostream>

#include <p8est_mesh.h>

namespace Gascoigne {

void
p8est_quadrant_init_fn(p8est_t* p8est,
                       p4est_topidx_t which_tree,
                       p8est_quadrant_t* quadrant)
{
  P8estMeshAgent::pforest_data_t* pforest_data =
    ((P8estMeshAgent::pforest_data_t*)(p8est->user_pointer));
  ((P8estMeshAgent::pquadrant_data_t*)(quadrant->p.user_data))->index =
    pforest_data->MAX_INDEX++;
}

template<>
P8estMeshAgent::~PForestMeshAgent()
{
  /* Destroy the p8est and the connectivity structure. */
  p8est_destroy(pforest);
  p8est_connectivity_destroy(conn);
}

template<>
void
P8estMeshAgent::read_inp(const std::string& gridname)
{
  conn = p8est_connectivity_read_inp(gridname.c_str());
  if (conn == NULL) {
    P4EST_LERRORF("Failed to read a valid connectivity from %s\n",
                  gridname.c_str());
    sc_abort();
  }

  /* Create a forest that is not refined; it consists of the root octant. */
  pforest = p8est_new(sc_MPI_COMM_NULL,
                      conn,
                      sizeof(pquadrant_data_t),
                      p8est_quadrant_init_fn,
                      &pforest_data);
}

template<>
void
P8estMeshAgent::write_vtk(const std::string& fname) const
{
  /* Write the forest to disk for visualization, one file per processor. */
  p8est_vtk_write_file(pforest, NULL, fname.c_str());
}

template<>
void
P8estMeshAgent::global_refine(IndexType refine_level)
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
  auto refine_fn = [](p8est_t* p8est,
                      p4est_topidx_t which_tree,
                      p8est_quadrant_t* quadrant) -> int { return 1; };
  /* Refine the forest iteratively, load balancing at each iteration.
   * This is important when starting with an unrefined forest */

  for (IndexType level = 0; level < refine_level; ++level) {
    p8est_refine(pforest, 0, refine_fn, p8est_quadrant_init_fn);
  }
}

template<>
void
P8estMeshAgent::refine_cells(IndexVector& ref)
{
  for (IndexType i = 0; i < pforest->trees->elem_count; ++i) {
    p8est_tree_t* tree = p8est_tree_array_index(pforest->trees, i);
    for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
      p8est_quadrant_t* quadrant =
        p8est_quadrant_array_index(&(tree->quadrants), j);
      pquadrant_data_t* data = ((pquadrant_data_t*)quadrant->p.user_data);
      if (std::find(ref.begin(), ref.end(), data->index) != ref.end()) {
        data->refine = true;
      } else {
        data->refine = false;
      }
    }
  }

  auto refine_fn = [](p8est_t* p8est,
                      p4est_topidx_t which_tree,
                      p8est_quadrant_t* quadrant) -> int {
    if (((pquadrant_data_t*)quadrant->p.user_data)->refine == true) {
      ((pquadrant_data_t*)quadrant->p.user_data)->refine = false;
      return 1;
    }
    return 0;
  };

  p8est_refine(pforest, 0, refine_fn, p8est_quadrant_init_fn);
}

template<>
IndexType
P8estMeshAgent::quad_count() const
{
  IndexType quad_count = 0;
  for (IndexType i = 0; i < trees_count(); ++i) {
    p8est_tree_t* tree = p8est_tree_array_index(pforest->trees, i);
    quad_count += tree->quadrants.elem_count;
  }
  return quad_count;
}

template class PForestMeshAgent<p8est,
                                p8est_tree_t,
                                p8est_quadrant_t,
                                p8est_connectivity_t>;

} // namespace Gascoigne