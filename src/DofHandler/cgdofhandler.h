/*----------------------------   cgdofhandler.h ---------------------------*/
/*      $Id:$                 */
#ifndef __cgdofhandler_H
#define __cgdofhandler_H
/*----------------------------   cgdofhandler.h ---------------------------*/

/**
 *
 * Copyright (C) 2019, 2020 by the Gascoigne 3D authors
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

#include "dofhandler.h"

/*-----------------------------------------*/

namespace Gascoigne {
/**
 * Handler for degrees of freedom for discretizations
 * with globally continuous functions, e.g. standard
 * lagrange discs.
 *
 * The DofHandler collects elements, each element has M^DIM
 * dofs, i.e. in 2D:
 *
 * M=2
 *
 * 2 - 3
 * |   |
 * 0 - 1
 *
 * M=3
 *
 * 6 - 7 - 8
 * |       |
 * 3   4   5
 * |       |
 * 0 - 1 - 2
 *
 * and general
 *
 * M(M-1) ..    M^2-1
 * .                .
 * .                .
 * M - ...      (2M-1)
 * 0 - 1 - ... - (M-1)
 *
 * In 3D, sorting in x-y-z direction
 *
 *   6----7
 *  /|   /|    y  z
 * 2----3 |    | /
 * | 4--|-5    |/
 * |/   |/     +---x
 * 0 -- 1
 *
 *
 * HANGING NODES
 *
 * 2:1 mesh refinement is allowed.
 *
 *
 ***/

/// DIM is spatial dimension, M dof's per element in each direction
template<int DIM, int M>
class CGDofHandler : public DofHandlerBase
{
protected:
  /// geometric coordinates for all dof's
  std::vector<Vertex<DIM>> nx;

public:
  CGDofHandler(){};
  ~CGDofHandler(){};

  /**
   * Initialize from mesh-structure in Base DofHandler<DIM> = GascoigneMesh<DIM>
   *
   * Each cell of the GascoigneMesh is filled with M^DIM dof's
   * dof's on nodes / edges / faces are shared
   **/
  void InitFromGascoigneMesh(const DofHandler<DIM>& GM);

  /// General Access
  IndexType dimension() const { return DIM; }
  IndexType dofs_per_element() const
  {
    return static_cast<IndexType>(pow(M, DIM));
  }
  IndexType ndofs() const { return nx.size(); }
  IndexType nelements() const { return nc.size() / dofs_per_element(); }
  IndexType nhanging() const { return HangingHandler.GetStructure()->size(); }

  IndexVector GetElement(IndexType iq) const
  {
    abort();

    IndexVector tmp;
    IndexType start = dofs_per_element() * iq;
    for (IndexType i = 0; i < dofs_per_element(); ++i)
      tmp.push_back(nc[start + i]);
    return tmp;
  }

  /// Access to coordinates
  std::vector<Vertex<DIM>>& GetVertexVector() { return nx; }
  const std::vector<Vertex<DIM>>& GetVertexVector() const { return nx; }

  const Vertex<DIM>& vertex(IndexType i) const { return nx[i]; }
  virtual const Vertex<2>& vertex2d(IndexType /*i*/) const { abort(); }
  virtual const Vertex<3>& vertex3d(IndexType /*i*/) const { abort(); }

  ////// Boundary
  const IndexVector* ElementOnBoundary(IndexType color) const
  {
    return &(BoundaryHandler.Cells(color));
  }
  const IndexVector* ElementLocalOnBoundary(IndexType color) const
  {
    return &(BoundaryHandler.Localind(color));
  }
  const IndexVector* ElementOnBoundary(IndexType /*degree*/,
                                       IndexType /*color*/) const
  {
    std::cerr << "Element on Boundary with degree not used!" << std::endl;
    abort();
  }
  const IndexVector* ElementLocalOnBoundary(IndexType /*degree*/,
                                            IndexType /*color*/) const
  {
    std::cerr << "ElementLocal on Boundary with degree not used!" << std::endl;
    abort();
  }

  IndexType VtkType(IndexType /*i*/) const { return (DIM == 2) ? 9 : 12; }

  IndexType vertex_of_cell(IndexType i, IndexType ii) const
  {
    return nc[dofs_per_element() * i + ii];
  }

  /// Dummy? Old Interface
  IndexType ncells() const { abort(); }
  IndexType nelements(IndexType /*degree*/) const { abort(); }
  IndexType nodes_per_cell(IndexType /*i*/) const { abort(); }
  IndexType nnodes() const { abort(); }

  IndexVector IndicesOfCell(IndexType /*iq*/) const
  {
    std::cerr << "CGDofHandler: Use GetElement" << std::endl;
    assert(0);
    return IndexVector();
  }

  IndexVector GetElement(IndexType /*degree*/, IndexType /*iq*/) const
  {
    std::cerr << "GetElement with degree not used" << std::endl;
    abort();
  }
};

template<>
inline const Vertex<2>&
CGDofHandler<2, 2>::vertex2d(IndexType i) const
{
  return vertex(i);
}
template<>
inline const Vertex<2>&
CGDofHandler<2, 3>::vertex2d(IndexType i) const
{
  return vertex(i);
}
template<>
inline const Vertex<2>&
CGDofHandler<2, 5>::vertex2d(IndexType i) const
{
  return vertex(i);
}
template<>
inline const Vertex<3>&
CGDofHandler<3, 2>::vertex3d(IndexType i) const
{
  return vertex(i);
}
template<>
inline const Vertex<3>&
CGDofHandler<3, 3>::vertex3d(IndexType i) const
{
  return vertex(i);
}
template<>
inline const Vertex<3>&
CGDofHandler<3, 5>::vertex3d(IndexType i) const
{
  return vertex(i);
}

} // namespace Gascoigne

/*----------------------------   dofhandler.h     ---------------------------*/
/* end of #ifndef __dofhandler_H */
#endif
/*----------------------------   dofhandler.h     ---------------------------*/
