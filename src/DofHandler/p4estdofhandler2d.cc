#ifndef P4_TO_P8
#include "p4estdofhandler2d.h"
#endif

#include <fstream>

#include "../Interface/gascoigne.h"

namespace Gascoigne {

/**
 * @brief Recreates the lnodes for current mesh refinement.
 *
 */
P4estDofHandler2d::P4estDofHandler2d(p4est_t* p4est, IndexType degree)
  : p4est(p4est)
{
  /* Create the ghost layer to learn about parallel neighbors. */
  p4est_ghost_t* ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
  /* Create a node numbering for continuous linear finite elements. */
  lnodes = p4est_lnodes_new(p4est, ghost, degree);
  /* Destroy the ghost structure -- no longer needed after node creation. */
  p4est_ghost_destroy(ghost);
  ghost = NULL;
}

P4estDofHandler2d::~P4estDofHandler2d()
{
  p4est_lnodes_destroy(lnodes);
}

/**
 * @param cell Index of the cell for witch the lnode index is looked up
 * @return IndexVector the lnode indices of the cell
 */
IndexVector
P4estDofHandler2d::get_nodes_of_cell(IndexType cell) const
{
  IndexVector nodes(P4EST_CHILDREN);
  for (IndexType i = 0; i < P4EST_CHILDREN; ++i) {
    nodes[i] = lnodes->element_nodes[P4EST_CHILDREN * cell + i];
  }
  return nodes;
}

/**
 * @return the number of lnodes in the mesh
 */
IndexType
P4estDofHandler2d::num_nodes() const
{
  return lnodes->num_local_elements;
}

#ifndef P4_TO_P8
/**
 * @brief
 *
 * @param file_name
 * @param vec
 */
void
P4estDofHandler2d::write_vtk(std::string file_name,
                             double time,
                             GhostVectorAgent& gva,
                             std::vector<std::string> vectors) const
{

  std::ofstream out(file_name.c_str());

  out << "# vtk DataFile Version 2.0 " << std::endl;
  out << "output from GascoigneStd" << std::endl;
  out << "ASCII" << std::endl;
  out << "DATASET UNSTRUCTURED_GRID" << std::endl;
  out << "FIELD FieldData 1" << std::endl;
  out << "TIME 1 1 double" << std::endl;
  out << time << std::endl;

  // Writing Points of the Mesh
  IndexType num_vertex = p4est->global_num_quadrants * 4;
  out << "POINTS " << num_vertex << " DOUBLE " << std::endl;

  for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);
    for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
      p4est_quadrant_t* quadrant =
        p4est_quadrant_array_index(&(tree->quadrants), j);
      for (IndexType k = 0; k < 4; ++k) {
        // Counting in a circle arount the Quad
        IndexType y = (k / 2) & 1;   // 0;0;1;1
        IndexType x = (y) ^ (k & 1); // 0;1;1;0

        double vxyz[3];
        double quad_lenght = P4EST_QUADRANT_LEN(quadrant->level);
        p4est_qcoord_to_vertex(p4est->connectivity,
                               i,
                               quadrant->x + x * quad_lenght,
                               quadrant->y + y * quad_lenght,
                               vxyz);
        Vertex2d coordinates(vxyz[0], vxyz[1], vxyz[2]);
        out << coordinates << " " << 0 << std::endl;
      }
    }
  }
  out << std::endl;

  // Writing mesh structur
  IndexType num_quads = p4est->global_num_quadrants;
  int lenght = num_quads * 5;

  out << std::endl << "CELLS " << num_quads << " " << lenght << std::endl;

  for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);
    for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
      IndexType id = (tree->quadrants_offset + j) * 4;
      int nle = 4;
      out << nle << " ";
      for (IndexType k = 0; k < nle; k++) {
        out << id + k << " ";
      }
      out << std::endl;
    }
  }

  out << std::endl << "CELL_TYPES " << num_quads << std::endl;
  for (int c = 0; c < num_quads; c++) {
    out << 9 << " ";
  }
  out << std::endl << std::endl;

  // Writing Vector
  out << "POINT_DATA " << num_vertex << std::endl;
  for (const std::string& vec_name : vectors) {
    GlobalVector* vec = gva[vec_name];
    if (!vec) {
      continue;
    }
    out << "SCALARS " << vec_name << " DOUBLE " << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (IndexType ind = 0; ind < num_quads; ind++) {
      for (IndexType j = 0; j < 4; ++j) {
        out << float((*vec)[lnodes->element_nodes[P4EST_CHILDREN * ind + j]])
            << std::endl;
      }
    }
    out << std::endl;
  }

  out.close();
}

#endif // P4_TO_P8

}