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

#include "p4estmeshagent.h"

#include <algorithm>
#include <iostream>

#include <p4est_mesh.h>

#include "../Common/dataformathandler.h"
#include "../Common/filescanner.h"

namespace Gascoigne {

P4estMeshAgent::P4estMeshAgent() {}
P4estMeshAgent::~P4estMeshAgent()
{
  /* Destroy the p4est and the connectivity structure. */
  p4est_destroy(p4est);
  p4est_connectivity_destroy(conn);
}

void
P4estMeshAgent::BasicInit(const ParamFile& pf)
{
  IndexType prerefine;
  std::string gridname;

  DataFormatHandler DFH;

  DFH.insert("gridname", &gridname, "none");
  DFH.insert("prerefine", &prerefine, 0);

  FileScanner FS(DFH);
  FS.readfile(pf, "Mesh");

  read_inp(gridname);
  global_refine(prerefine);
}

void
init_fn(p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* quadrant)
{
  P4estMeshAgent::p4est_data_t* p4est_data =
    ((P4estMeshAgent::p4est_data_t*)(p4est->user_pointer));
  ((P4estMeshAgent::quadrant_data_t*)(quadrant->p.user_data))->index =
    p4est_data->MAX_INDEX++;
}

void
P4estMeshAgent::read_inp(const std::string& gridname)
{
  conn = p4est_connectivity_read_inp(gridname.c_str());
  if (conn == NULL) {
    P4EST_LERRORF("Failed to read a valid connectivity from %s\n",
                  gridname.c_str());
    sc_abort();
  }

  /* Create a forest that is not refined; it consists of the root octant. */
  p4est = p4est_new(
    sc_MPI_COMM_NULL, conn, sizeof(quadrant_data_t), init_fn, &p4est_data);
}

void
P4estMeshAgent::write_vtk(const std::string& fname) const
{
  /* Write the forest to disk for visualization, one file per processor. */
  p4est_vtk_write_file(p4est, NULL, fname.c_str());
}

void
P4estMeshAgent::global_refine(IndexType refine_level)
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
    p4est_refine(p4est, 0, refine_fn, init_fn);
  }
}

void
P4estMeshAgent::refine_cells(IndexVector& ref)
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
  for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);
    for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
      p4est_quadrant_t* quadrant =
        p4est_quadrant_array_index(&(tree->quadrants), j);
      quadrant_data_t* data = ((quadrant_data_t*)quadrant->p.user_data);
      if (std::find(ref.begin(), ref.end(), data->index) != ref.end()) {
        data->refine = true;
      } else {
        data->refine = false;
      }
    }
  }

  auto refine_fn = [](p4est_t* p4est,
                      p4est_topidx_t which_tree,
                      p4est_quadrant_t* quadrant) -> int {
    if (((quadrant_data_t*)quadrant->p.user_data)->refine == true) {
      ((quadrant_data_t*)quadrant->p.user_data)->refine = false;
      return 1;
    }
    return 0;
  };

  p4est_refine(p4est, 0, refine_fn, init_fn);
}

} // namespace Gascoigne