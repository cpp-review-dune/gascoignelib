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

// namespace Gascoigne {
// /**
//  * @brief
//  *
//  * @param file_name
//  * @param vec
//  */
// void
// P4estDofHandler3d::write_vtk(std::string file_name,
//                              double time,
//                              GhostVectorAgent& gva,
//                              std::vector<std::string> vectors) const
// {

//   std::ofstream out(file_name.c_str());

//   out << "# vtk DataFile Version 2.0 " << std::endl;
//   out << "output from GascoigneStd" << std::endl;
//   out << "ASCII" << std::endl;
//   out << "DATASET UNSTRUCTURED_GRID" << std::endl;
//   out << "FIELD FieldData 1" << std::endl;
//   out << "TIME 1 1 double" << std::endl;
//   out << time << std::endl << std::endl;

//   // Writing Points of the Mesh
//   IndexType num_vertex = p4est->global_num_quadrants * P8EST_CHILDREN;
//   out << "POINTS " << num_vertex << " DOUBLE " << std::endl;

//   for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
//     p8est_tree_t* tree = p8est_tree_array_index(p4est->trees, i);
//     for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
//       p8est_quadrant_t* quadrant =
//         p8est_quadrant_array_index(&(tree->quadrants), j);
//       for (IndexType k = 0; k < 8; ++k) {
//         // Counting in a circle arount the Quad
//         IndexType y = (k / 2) % 2; // 0;0;1;1;0;0;1;1
//         IndexType x = k % 2;       // 0;1;1;0;0;1;1;0
//         IndexType z = k / 4;       // 0;0;0;0;1;1;1;1

//         double vxyz[3];
//         double quad_lenght = P8EST_QUADRANT_LEN(quadrant->level);
//         p8est_qcoord_to_vertex(p4est->connectivity,
//                                i,
//                                quadrant->x + x * quad_lenght,
//                                quadrant->y + y * quad_lenght,
//                                quadrant->z + z * quad_lenght,
//                                vxyz);
//         Vertex3d coordinates(vxyz[0], vxyz[1], vxyz[2]);
//         out << coordinates << std::endl;
//       }
//     }
//   }
//   out << std::endl;

//   // Writing mesh structur
//   IndexType num_quads = p4est->global_num_quadrants;
//   int lenght = num_quads * (P8EST_CHILDREN + 1);

//   out << std::endl << "CELLS " << num_quads << " " << lenght << std::endl;

//   for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
//     p8est_tree_t* tree = p8est_tree_array_index(p4est->trees, i);
//     for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
//       IndexType id = (tree->quadrants_offset + j) * P8EST_CHILDREN;
//       int nle = P8EST_CHILDREN;
//       out << nle << " ";
//       for (IndexType k = 0; k < nle; k++) {
//         out << id + k << " ";
//       }
//       out << std::endl;
//     }
//   }

//   out << std::endl << "CELL_TYPES " << num_quads << std::endl;
//   for (int c = 0; c < num_quads; c++) {
//     out << 11 << " ";
//   }
//   out << std::endl << std::endl;

//   // Writing Vector
//   out << "POINT_DATA " << num_vertex << std::endl;
//   for (const std::string& vec_name : vectors) {
//     GlobalVector* vec = gva[vec_name];
//     if (!vec) {
//       continue;
//     }
//     out << "SCALARS " << vec_name << " DOUBLE " << std::endl;
//     out << "LOOKUP_TABLE default" << std::endl;
//     for (IndexType ind = 0; ind < num_quads; ind++) {
//       for (IndexType j = 0; j < P8EST_CHILDREN; ++j) {
//         out << float((*vec)[lnodes->element_nodes[P8EST_CHILDREN * ind + j]])
//             << std::endl;
//       }
//     }
//     out << std::endl;
//   }

//   out.close();
// }

// }