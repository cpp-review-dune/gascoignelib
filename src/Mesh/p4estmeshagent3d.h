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

#ifndef __p4estmesh3d_h
#define __p4estmesh3d_h

#include <set>
#include <string>
#include <utility>
#include <vector>

#include <p8est.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_vtk.h>

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
#include "p4estmeshagentbase.h"

/*---------------------------------------------------*/

namespace Gascoigne {

class P4estMeshAgent3d : public P4estMeshAgentBase
{
private:
  p8est_t* pforest;
  p8est_connectivity_t* conn;
  p8est_lnodes_t* plnodes;

  pforest_data_t pforest_data;

public:
  P4estMeshAgent3d(const std::string& gridname,
                   IndexType prerefine = 0,
                   IndexType comp = 1);
  virtual ~P4estMeshAgent3d();

  virtual IndexType trees_count() const;
  virtual IndexType quad_count() const;

  virtual void write_vtk(const std::string& fname) const;
  virtual void global_refine(IndexType n = 1);
  virtual void refine_cells(IndexVector& ref);
};

} // namespace Gascoigne

/*---------------------------------------------------*/

#endif //__p4estmesh3d_h
