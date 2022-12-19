/*----------------------------   p4estdofhandler.h ---------------------------*/
/*      $Id:$                 */
#ifndef __p4estdofhandler3d_H
#define __p4estdofhandler3d_H
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

#include <p8est_lnodes.h>

#include "../Interface/gascoigne.h"

#include "p4estdofhandler.h"

namespace Gascoigne {

class P4estDofHandler3d : public P4estDofHandler
{
private:
  p8est_lnodes_t* lnodes;

public:
  P4estDofHandler3d(p8est_t* pforest, IndexType degree);
  virtual ~P4estDofHandler3d();

  IndexVector get_nodes_of_cell(IndexType cell) const;
  IndexType num_nodes() const;
};

}
#endif //__p4estdofhandler3d_H