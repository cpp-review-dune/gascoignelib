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
  
  reorder_hanging_nodes();
  generate_lnode_pos();
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
    nodes[i] = get_node_of_cell( cell, i);
  }
  return nodes;
}

/**
 * @param cell Index of the cell for witch the lnode index is looked up
 * @return IndexVector the lnode indices of the cell
 */
IndexType
P4estDofHandler2d::get_node_of_cell(IndexType cell, IndexType i) const
{
  return lnodes->element_nodes[nodes_per_cell() * cell + i];
}

/**
 * @return the number of lnodes in the mesh
 */
IndexType
P4estDofHandler2d::num_nodes() const
{
  return lnodes->num_local_nodes + hn.size();
}

/**
 * @return the number of lnodes in the mesh
 */
IndexType
P4estDofHandler2d::num_haning() const
{
  return hn.size();
}


std::array<MatrixEntryType, DIM>
P4estDofHandler2d::vertex(IndexType node_index) const {
  return lnode_pos[node_index];
}


// TODO: Impelement for Degree = 2
static const int    ones = P4EST_CHILDREN - 1;  /**< One bit per dimension. */
/** Decode the information from p{4,8}est_lnodes_t for a given element.
 *
 * \see p4est_lnodes.h for an in-depth discussion of the encoding.
 * \param [in] face_code         Bit code as defined in p{4,8}est_lnodes.h.
 * \param [out] hanging_corner   Undefined if no node is hanging.
 *                               If any node is hanging, this contains
 *                               one integer per corner, which is -1
 *                               for corners that are not hanging,
 *                               and the number of the non-hanging
 *                               corner on the hanging face/edge otherwise.
 *                               For faces in 3D, it is diagonally opposite.
 * \return true if any node is hanging, false otherwise.
 */
static int
lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN])
{
  if (face_code) {
    const int           c = (int) (face_code & ones);
    int                 i, h;
    int                 work = (int) (face_code >> P4EST_DIM);

    /* These two corners are never hanging by construction. */
    hanging_corner[c] = hanging_corner[c ^ ones] = -1;
    for (i = 0; i < P4EST_DIM; ++i) {
      /* Process face hanging corners. */
      h = c ^ (1 << i);
      hanging_corner[h ^ ones] = (work & 1) ? c : -1;
#ifdef P4_TO_P8
      /* Process edge hanging corners. */
      hanging_corner[h] = (work & P4EST_CHILDREN) ? c : -1;
#endif
      work >>= 1;
    }
    return 1;
  }
  return 0;
}

const IndexType corner_degree1to2[8] = {0,2,6,8,18,20,24,26};
void P4estDofHandler2d::insert_hn(IndexType global_quad_id, IndexType c1, IndexType c2, IndexType where){
  if(_degree == 2){
    c1 = corner_degree1to2[c1];
    c2 = corner_degree1to2[c2];
  }

  IndexType n0 = get_node_of_cell(global_quad_id, c1);
  IndexType n1 = get_node_of_cell(global_quad_id, c2);

  // Determen if hn exists or add hn and get id
  IndexType hn_id = 0;
  for(IndexType j = 0; j < hn.size(); ++j){
    if((hn[j][0] == n0 && hn[j][1] == n1) || (hn[j][0] == n1 && hn[j][1] == n0)) {
      hn_id = lnodes->num_local_nodes + j;
      break;
    }
  }
  if(!hn_id){
    hn_id = num_nodes();
    std::vector<IndexType> node = {n0, n1};
    hn.emplace_back(node);
  }

  // Adding the Hanging Node Indices to lnode strukture !! Do NOT do this at home !!
  lnodes->element_nodes[nodes_per_cell() * global_quad_id + where] = hn_id;
}

void
P4estDofHandler2d::reorder_hanging_nodes(){
  for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);
    for (IndexType q = 0; q < tree->quadrants.elem_count; ++q) {
      IndexType global_quad_id = (tree->quadrants_offset + q);
      int hanging_corner[P4EST_CHILDREN];
      if(!lnodes_decode2 (lnodes->face_code[global_quad_id], hanging_corner)){
        continue;
      }

      if(_degree == 2){
        for(IndexType k = 0; k < P4EST_CHILDREN; ++k){
          if(hanging_corner[k] != -1){
            IndexType between = (corner_degree1to2[k] + corner_degree1to2[hanging_corner[k]]) / 2;
            lnodes->element_nodes[nodes_per_cell() * global_quad_id + corner_degree1to2[k]] = get_node_of_cell(global_quad_id, between);
            insert_hn(global_quad_id, k, hanging_corner[k], between);
          }
        }

        for(IndexType c1 = 0; c1 < P4EST_CHILDREN-1; ++c1){
          if(hanging_corner[c1] != -1){
            for(IndexType c2 = c1+1; c2 < P4EST_CHILDREN; ++c2){
              if(hanging_corner[c2] != -1){
                IndexType between = (corner_degree1to2[c1] + corner_degree1to2[c2]) / 2;
                insert_hn(global_quad_id, c1, c2, between);
              }
            }
          }
        }
      } else {
        for(IndexType k = 0; k < P4EST_CHILDREN; ++k){
          if(hanging_corner[k] != -1){
            insert_hn(global_quad_id, k, hanging_corner[k], k);
          }
        }
        
      }
    }
  }
}

void
P4estDofHandler2d::generate_lnode_pos()
{
  lnode_pos.resize(num_nodes());

  for (IndexType i = 0; i < p4est->trees->elem_count; ++i) {
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);
    for (IndexType q = 0; q < tree->quadrants.elem_count; ++q) {
      IndexType global_quad_id = (tree->quadrants_offset + q);
      p4est_quadrant_t* quadrant =
        p4est_quadrant_array_index(&(tree->quadrants), q);

      double _degree_length = P4EST_QUADRANT_LEN(quadrant->level) / _degree; //< Distance between two lnodes at _degree
      IndexType degp1 = _degree + 1; //< Degree Plus One for Indexing

      for (IndexType k = 0; k < nodes_per_cell(); ++k) {
        double vxyz[3] = { 0, 0, 0 };
        p4est_qcoord_to_vertex(p4est->connectivity,
                               i,
                               quadrant->x + k % (degp1) * _degree_length,
                               quadrant->y + (k / degp1) % degp1 * _degree_length,
#ifdef P4_TO_P8
                               quadrant->z + (k / (degp1 * degp1)) % degp1 * _degree_length,
#endif
                               vxyz);

        lnode_pos[get_node_of_cell( global_quad_id, k)] = {
          vxyz[0], vxyz[1],
#ifdef P4_TO_P8
           vxyz[2]
#endif
        };
      }
    }
  }
  for(IndexType i = lnodes->num_local_nodes; i < num_nodes(); ++i){
    auto node = hn[i - lnodes->num_local_nodes];
    lnode_pos[i] = {
      (lnode_pos[node[0]][0] + lnode_pos[node[1]][0]) * 0.5,
      (lnode_pos[node[0]][1] + lnode_pos[node[1]][1]) * 0.5,
#ifdef P4_TO_P8
      (lnode_pos[node[0]][2] + lnode_pos[node[1]][2]) * 0.5
#endif
    };
  }
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
  out << "POINTS " << num_nodes() << " DOUBLE " << std::endl;

  for (IndexType i = 0; i < num_nodes(); ++i) {
    out << lnode_pos[i][0] << " " << lnode_pos[i][1] << " ";
#ifdef P4_TO_P8
    out << lnode_pos[i][2];
#else
    out << "0";
#endif
    out << std::endl;
  }
  out << std::endl;

  // std::cout << "# Hanging Nodes\n";
  // for(IndexType i = 0; i < hn.size(); ++i){
  //   std::cout << "# " << (i + lnodes->num_local_nodes) << " " << hn[i][0] << " " << hn[i][1] << std::endl;
  // }

  // Writing mesh structur
  IndexType num_cells = p4est->global_num_quadrants;
  int lenght = num_cells * (nodes_per_cell() + 1);

  out << std::endl << "CELLS " << num_cells << " " << lenght << std::endl;

  IndexType* node_order = nullptr;
#ifdef P4_TO_P8
  IndexType node_order_q13[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
  IndexType node_order_q23[27] = { 0,  2,  8,  6,  18, 20, 26, 24, 1,
                                   5,  7,  3,  19, 23, 25, 21, 9,  11,
                                   17, 15, 12, 14, 10, 16, 4,  22, 13 };
  if (_degree == 1) {
    node_order = node_order_q13;
  } else if (_degree == 2) {
    node_order = node_order_q23;
  }
#else
  IndexType node_order_q12[4] = { 0, 1, 3, 2 };
  IndexType node_order_q22[9] = { 0, 2, 8, 6, 1, 5, 7, 3, 4 };
  if (_degree == 1) {
    node_order = node_order_q12;
  } else if (_degree == 2) {
    node_order = node_order_q22;
  }
#endif

  for (IndexType ind = 0; ind < num_cells; ind++) {
    int nle = nodes_per_cell();
    out << nle << " ";
    for (IndexType j = 0; j < nodes_per_cell(); ++j) {
      out << get_node_of_cell( ind, node_order[j])
          << " ";
    }
    out << std::endl;
  }

  out << std::endl << "CELL_TYPES " << num_cells << std::endl;
  IndexType celltype = -1;
#ifdef P4_TO_P8
  if (_degree == 1) {
    celltype = 12;
  } else if (_degree == 2) {
    celltype = 29;
  }
#else
  if (_degree == 1) {
    celltype = 9;
  } else if (_degree == 2) {
    celltype = 28;
  }
#endif
  for (int c = 0; c < num_cells; c++) {
    out << celltype << " ";
  }
  out << std::endl << std::endl;

  // Writing Vector
  out << "POINT_DATA " << num_nodes() << std::endl;
  for (const std::string& vec_name : vectors) {
    GlobalVector* vec = gva[vec_name];
    if (!vec) {
      continue;
    }
    out << "SCALARS " << vec_name << " DOUBLE " << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (IndexType ind = 0; ind < num_nodes(); ind++) {
      out << float((*vec)[ind]) << std::endl;
    }
    out << std::endl;
  }

  out.close();
}

}