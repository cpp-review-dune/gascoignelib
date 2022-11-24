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

#include "p4estmesh2d.h"

#include <p4est_mesh.h>

namespace Gascoigne {

void
P4estMesh2d::read_inp(const std::string& gridname)
{
  conn = p4est_connectivity_read_inp(gridname.c_str());
  if (conn == NULL) {
    P4EST_LERRORF("Failed to read a valid connectivity from %s\n",
                  gridname.c_str());
    sc_abort();
  }
  /* Create a forest that is not refined; it consists of the root octant. */
  p4est = p4est_new(sc_MPI_COMM_NULL, conn, 0, NULL, NULL);
}

void
P4estMesh2d::refine(const IntVector& cell_ref,
                    const IntVector& cell_coarse){ TO_DO }

std::set<int> P4estMesh2d::CellNeighbours(int iq) const
{
  std::set<int> neighbors;
  p4est_mesh_face_neighbor_t iterator;

  p4est_mesh_face_neighbor_init2(&iterator, p4est, NULL, NULL, 0, iq);
  TO_DO
  return std::set<int>();
}

} // namespace Gascoigne