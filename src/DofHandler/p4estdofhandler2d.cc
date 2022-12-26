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
  : P4estDofHandler(P4EST_DIM, degree)
  , p4est(p4est)
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
  IndexVector nodes(nodes_per_cell());
  for (IndexType i = 0; i < nodes_per_cell(); ++i) {
    nodes[i] = lnodes->element_nodes[nodes_per_cell() * cell + i];
  }
  return nodes;
}

/**
 * @return the number of lnodes in the mesh
 */
IndexType
P4estDofHandler2d::num_nodes() const
{
  return lnodes->num_local_nodes;
}

/**
 * @brief Writes a vtk file for further processing
 *
 * @param file_name Where to warite
 * @param time What timestep to date it
 * @param gva A GhostVectorAgent from the dof
 * @param vectors Vectors from the gva that are written
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
  out << time << std::endl << std::endl;

  // Writing Points of the Mesh
  IndexType num_vertex = p4est->global_num_quadrants * nodes_per_cell();
  out << "POINTS " << num_vertex << " DOUBLE " << std::endl;

  for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);
    for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
      p4est_quadrant_t* quadrant =
        p4est_quadrant_array_index(&(tree->quadrants), j);
      for (IndexType k = 0; k < nodes_per_cell(); ++k) {
        // Counting in a circle arount the Quad
#ifdef P4_TO_P8
        IndexType z = (k / 4) % 2; // 0;0;0;0;1;1;1;1
#endif
        IndexType y = (k / 2) % 2; // 0;0;1;1;0;0;1;1
        IndexType x = k % 2;       // 0;1;1;0;0;1;1;0

        double vxyz[3] = { 0, 0, 0 };
        double quad_lenght = P4EST_QUADRANT_LEN(quadrant->level);
        p4est_qcoord_to_vertex(p4est->connectivity,
                               i,
                               quadrant->x + x * quad_lenght,
                               quadrant->y + y * quad_lenght,
#ifdef P4_TO_P8
                               quadrant->z + z * quad_lenght,
#endif
                               vxyz);
        Vertex3d coordinates(vxyz[0], vxyz[1], vxyz[2]);
        out << coordinates << std::endl;
      }
    }
  }
  out << std::endl;

  // Writing mesh structur
  IndexType num_cells = p4est->global_num_quadrants;
  int lenght = num_cells * (P4EST_CHILDREN + 1);

  out << std::endl << "CELLS " << num_cells << " " << lenght << std::endl;

  for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);
    for (IndexType j = 0; j < tree->quadrants.elem_count; ++j) {
      IndexType id = (tree->quadrants_offset + j) * P4EST_CHILDREN;
      int nle = P4EST_CHILDREN;
      out << nle << " ";
      for (IndexType k = 0; k < nle; k++) {
        out << id + k << " ";
      }
      out << std::endl;
    }
  }

  out << std::endl << "CELL_TYPES " << num_cells << std::endl;
#ifndef P4_TO_P8
  IndexType celltype = 8;
#else
  IndexType celltype = 11;
#endif
  for (int c = 0; c < num_cells; c++) {
    out << celltype << " ";
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
    for (IndexType ind = 0; ind < num_cells; ind++) {
      for (IndexType j = 0; j < nodes_per_cell(); ++j) {
        out << float((*vec)[lnodes->element_nodes[nodes_per_cell() * ind + j]])
            << std::endl;
      }
    }
    out << std::endl;
  }

  out.close();
}

}