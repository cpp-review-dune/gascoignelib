#ifndef P4_TO_P8
#include "p4estdofhandler2d.h"
#endif

#include "../Interface/gascoigne.h"

namespace Gascoigne {

/**
 * @brief Recreates the lnodes for current mesh refinement.
 *
 */
P4estDofHandler2d::P4estDofHandler2d(p4est_t* p4est, IndexType degree)
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

}