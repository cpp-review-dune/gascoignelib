/*----------------------------   p4estdofhandler.h ---------------------------*/
/*      $Id:$                 */
#ifndef __p4estdofhandlerbase_H
#define __p4estdofhandlerbase_H
/*----------------------------   p4estdofhandler.h ---------------------------*/

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

#include "../Common/compvector.h"
#include "../Interface/gascoigne.h"
#include "../LinAlg/columnstencil.h"
#include "../Solver/ghostagent.h"

namespace Gascoigne {

class P4estDofHandlerBase
{
protected:
  ColumnStencil _hn;

  P4estDofHandlerBase(){};

public:
  IndexType nodes_per_cell() const { return pow(degree() + 1, dimension()); };
  virtual IndexVector get_nodes_of_cell(IndexType cell) const = 0;
  virtual IndexType get_node_of_cell(IndexType cell, IndexType i) const = 0;

  virtual IndexType nnodes() const = 0;
  virtual IndexType nnothanging() const = 0;
  virtual IndexType nhanging() const = 0;
  IndexType is_haning(IndexType index) const { return index > nnothanging(); };

  virtual IndexType dimension() const = 0;
  virtual IndexType degree() const = 0; //< degree of the lnodes

  virtual void write_vtk(std::string file_name,
                         double time,
                         GhostVectorAgent& gva,
                         std::vector<std::string> vectors) const = 0;
};

}
#endif //__p4estdofhandlerbase_H