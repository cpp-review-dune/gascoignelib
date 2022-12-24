/*----------------------------   p4estdofhandler.h ---------------------------*/
/*      $Id:$                 */
#ifndef __p4estdofhandler2d_H
#define __p4estdofhandler2d_H
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

#include <p4est.h>
#include <p4est_lnodes.h>

#include "../Common/vertex.h"
#include "../Interface/gascoigne.h"

#include "p4estdofhandler.h"

namespace Gascoigne {

class P4estDofHandler2d : public P4estDofHandler
{
private:
  p4est_t* p4est;
  p4est_lnodes_t* lnodes;

public:
  P4estDofHandler2d(p4est_t* pforest, IndexType degree);
  virtual ~P4estDofHandler2d();

  IndexVector get_nodes_of_cell(IndexType cell) const;
  IndexType num_nodes() const;
  IndexType dimension() const { return 2; };
  Vertex2d vertex2d(IndexType node_index) const;

  void write_vtk(std::string file_name, GlobalVector vec) const;
};
}
#endif //__p4estdofhandler2d_H