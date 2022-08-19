/**
 *
 * Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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

#ifndef __ColumnStencil_h
#define __ColumnStencil_h

#include "gascoigne.h"
#include "sparsestructureinterface.h"
#include "stencilinterface.h"

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  ... comments ColumnStencil

////
////
/////////////////////////////////////////////

class ColumnStencil : public virtual StencilInterface
{
protected:
  IndexVector scol, sstart;

public:
  //
  ////  Con(De)structor
  //

  ColumnStencil()
    : StencilInterface()
  {}
  ~ColumnStencil() {}

  const IndexVector& col() const { return scol; }
  IndexVector& col() { return scol; }
  const IndexVector& start() const { return sstart; }
  IndexVector& start() { return sstart; }

  IndexType n() const { return sstart.size() - 1; }
  IndexType nentries() const { return scol.size(); }
  IndexType rowsize(IndexType i) const
  {
    assert(i + 1 < sstart.size());
    return sstart[i + 1] - sstart[i];
  }

  IndexType& col(IndexType pos)
  {
    assert(pos < scol.size());
    return scol[pos];
  }
  const IndexType& col(IndexType pos) const
  {
    assert(pos < scol.size());
    return scol[pos];
  }
  IndexType& start(IndexType i)
  {
    assert(i < sstart.size());
    return sstart[i];
  }
  const IndexType& start(IndexType i) const
  {
    assert(i < sstart.size());
    return sstart[i];
  }
  IndexType& stop(IndexType i)
  {
    assert(i + 1 < sstart.size());
    return sstart[i + 1];
  }
  const IndexType& stop(IndexType i) const
  {
    assert(i + 1 < sstart.size());
    return sstart[i + 1];
  }

  void memory(const SparseStructureInterface*);
  void memory(int n, int nt);

  virtual int Find(int i, int j) const
  {
    for (int pos = start(i); pos < stop(i); pos++) {
      if (col(pos) == j)
        return pos;
    }
    std::cerr << "UnstructuredStencil::Find()";
    std::cerr << "no such coupling: " << i << " " << j << std::endl;
    abort();
    return -1;
  }

  std::ostream& Write(std::ostream& os) const;
};
} // namespace Gascoigne

#endif
