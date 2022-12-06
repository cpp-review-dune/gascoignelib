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
#include <p4est_vtk.h>
#include <p8est.h>
#include <p8est_vtk.h>

#include "../Common/compvector.h"
#include "../Common/paramfile.h"
#include "../Common/triple.h"
#include "../Common/vertex.h"
#include "../Interface/gascoigne.h"

#include "curvedshapes.h"
#include "edge.h"
#include "hanglist.h"
#include "hierarchicalmesh2d.h"

/*---------------------------------------------------*/

namespace Gascoigne {

template<typename PForest, typename PTree, typename PQuad, typename PConn>
class PForestMeshAgent
{
public:
  struct pquadrant_data_t
  {
    IndexType index;
    bool refine; ///< flag when true gets refinde in refine_cells
  };

  struct pforest_data_t
  {
    IndexType MAX_INDEX = 0;
  };

private:
  PForest* pforest;
  PConn* conn;
  pforest_data_t pforest_data;

public:
  PForestMeshAgent();
  virtual ~PForestMeshAgent();

  virtual void BasicInit(const ParamFile& pf);

  virtual IndexType trees_count() const;
  virtual IndexType quad_count() const;

  virtual void read_inp(const std::string& fname);
  virtual void write_vtk(const std::string& fname) const;
  virtual void global_refine(IndexType n);
  virtual void refine_cells(IndexVector& ref);
};

template<typename PForest, typename PTree, typename PQuad, typename PConn>
inline IndexType
PForestMeshAgent<PForest, PTree, PQuad, PConn>::trees_count() const
{
  return pforest->trees->elem_count;
}

using P4estMeshAgent =
  PForestMeshAgent<p4est, p4est_tree_t, p4est_quadrant_t, p4est_connectivity_t>;
using P8estMeshAgent =
  PForestMeshAgent<p8est, p8est_tree_t, p8est_quadrant_t, p8est_connectivity_t>;

} // namespace Gascoigne

/*---------------------------------------------------*/

#endif //__p4estmesh2d_h
