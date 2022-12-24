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

#include "p4estdofhandler3d.h"

#include <p4est_to_p8est.h>

#define P4estDofHandler2d P4estDofHandler3d

#include "p4estdofhandler2d.cc"

#undef P4estMeshAgent

namespace Gascoigne {
/**
 * @brief
 *
 * @param file_name
 * @param vec
 */
void
P4estDofHandler2d::write_vtk(std::string file_name, GlobalVector vec) const
{
  IndexType time = 0;
  std::ofstream out(file_name.c_str());

  out << "# vtk DataFile Version 2.0 " << std::endl;
  out << "output from GascoigneStd" << std::endl;
  out << "ASCII" << std::endl;
  out << "DATASET UNSTRUCTURED_GRID" << std::endl;
  out << "FIELD FieldData 1" << std::endl;
  out << "TIME 1 1 double" << std::endl;
  out << time << std::endl;

  // Writing Points of the Mesh
  IndexType nn = lnodes->num_local_elements;
  out << "POINTS " << nn << " DOUBLE " << std::endl;

  for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
    p8est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);
    for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
      p8est_quadrant_t* quadrant =
        p4est_quadrant_array_index(&(tree->quadrants), j);
      double vxyz[3];
      p8est_qcoord_to_vertex(
        p4est->connectivity, i, quadrant->x, quadrant->y, quadrant->z, vxyz);
      Vertex2d coordinates(vxyz[0], vxyz[1], vxyz[2]);
      out << coordinates << " " << 0 << std::endl;
    }
  }
  out << std::endl;

  // Writing mesh structur
//   IndexType ne = p4est->global_num_quadrants;
//   int lenght = ne*5;

//   out << std::endl << "CELLS " << ne << " " << lenght << std::endl;

//   for (int c = 0; c < ne; c++) {
//     int nle = 4;
//     out << nle << " ";
//     for (int ii = 0; ii < nle; ii++) {
//       out << mesh->vertex_of_cell(c, ii) << " ";
//     }
//     out << endl;
//   }
//   out << endl << "CELL_TYPES " << ne << endl;
//   for (int c = 0; c < ne; c++) {
//     out << mesh->VtkType(c) << " ";
//   }
//   out << endl;

  out.close();
}

}