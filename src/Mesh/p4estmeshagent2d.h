/**
 *
 * Copyright (C) 2004, 2005, 2006, 2011 by the Gascoigne 3D authors
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

#ifndef __p4estmesh2d_h
#define __p4estmesh2d_h

#include <set>
#include <string>
#include <utility>
#include <vector>

#include <p4est.h>
#include <p4est_connectivity.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>

#include "../Common/compvector.h"
#include "../Common/dataformathandler.h"
#include "../Common/filescanner.h"
#include "../Common/paramfile.h"
#include "../Common/triple.h"
#include "../Common/vertex.h"
#include "../Interface/gascoigne.h"

#include "curvedshapes.h"
#include "edge.h"
#include "hanglist.h"
#include "hierarchicalmesh2d.h"
#include "p4estmeshagent.h"

/*---------------------------------------------------*/

namespace Gascoigne {

class P4estMeshAgent2d : public P4estMeshAgent
{
private:
  p4est_t* pforest;
  p4est_connectivity_t* conn;

public:
  P4estMeshAgent2d(const std::string& gridname,
                   IndexType prerefine = 0,
                   IndexType comp = 1);
  virtual ~P4estMeshAgent2d();

  IndexType num_cells() const;

  void write_vtk(const std::string& fname) const;
  void global_refine(IndexType n = 1);
  void refine_cells(IndexVector& ref);

  std::shared_ptr<P4estDofHandlerBase> create_dofhandler(
    IndexType degree) const;
};

} // namespace Gascoigne

/*---------------------------------------------------*/

#endif //__p4estmesh2d_h
