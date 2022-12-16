/*----------------------------   dofhandler.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dofhandler_H
#define __dofhandler_H
/*----------------------------   dofhandler.h     ---------------------------*/

/**
 *
 * Copyright (C) 2018 by the Gascoigne 3D authors
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

#include "../Common/paramfile.h"
#include "../Common/vertex.h"
#include "../Interface/gascoigne.h"
#include "../Mesh/boundaryindexhandler.h"
#include "../Mesh/hangingindexhandler.h"
#include "../Mesh/patchindexhandler.h"

#include "dofhandlerbase.h"

/*-----------------------------------------*/

namespace Gascoigne {
template<int DIM>
class DofHandler : public DofHandlerBase
{
protected:
  // basic
  std::vector<Vertex<DIM>> nx;

public:
  DofHandler(){};
  ~DofHandler(){};

  std::string GetName() const;
  std::vector<Vertex<DIM>>& GetVertexVector() { return nx; }
  const std::vector<Vertex<DIM>>& GetVertexVector() const { return nx; }

  IndexType dimension() const { return DIM; }
  IndexType nnodes() const { return nx.size(); }
  IndexType ncells() const { return nc.size() / ((DIM == 2) ? 4 : 8); }
  IndexType nelements(IndexType degree) const
  {
    if (degree == 1)
      return ncells();
    else if (degree == 2)
      return npatches();
    else
      abort();
  }

  IndexType nhanging() const { return HangingHandler.GetStructure()->size(); }

  IndexType nodes_per_cell(IndexType /*i*/) const { return (DIM == 2) ? 4 : 8; }
  IndexType nodes_per_element(IndexType degree) const
  {
    if (degree == 1)
      return nodes_per_cell(0);
    else if (degree == 2)
      return nodes_per_patch();
    else
      abort();
  }

  IndexType VtkType(IndexType /*i*/) const { return (DIM == 2) ? 9 : 12; }

  const Vertex<DIM>& vertex(IndexType i) const
  {
    assert(i < nx.size());
    return nx[i];
  }

  virtual const Vertex<2>& vertex2d(IndexType i) const;
  virtual const Vertex<3>& vertex3d(IndexType i) const;

  IndexType vertex_of_cell(IndexType i, IndexType ii) const
  {
    assert(nodes_per_cell(i) * i + ii < nc.size());
    return nc[nodes_per_cell(i) * i + ii];
  }

  virtual IndexType CornerIndices(IndexType degree,
                                  IndexType iq,
                                  IndexType ii) const
  {
    if (degree == 1) {
      assert(nodes_per_cell(iq) * iq + ii < nc.size());
      return nc[nodes_per_cell(iq) * iq + ii];
    } else if (degree == 2) {
      return PatchHandler.CoarseIndices(iq)[ii];
    }
    assert(0);
    return 0;
  }

  IndexVector IndicesOfCell(IndexType iq) const;

  IndexVector GetElement(IndexType degree, IndexType iq) const
  {
    if (degree == 1)
      return IndicesOfCell(iq);
    else if (degree == 2)
      return *IndicesOfPatch(iq);
    abort();
  }
};

} // namespace Gascoigne

/*----------------------------   dofhandler.h     ---------------------------*/
/* end of #ifndef __dofhandler_H */
#endif
/*----------------------------   dofhandler.h     ---------------------------*/
