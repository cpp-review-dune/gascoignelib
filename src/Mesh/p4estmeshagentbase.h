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

#ifndef __p4estmeshagentbase_h
#define __p4estmeshagentbase_h

#include <memory>
#include <string>

#include "../Common/dataformathandler.h"
#include "../Common/filescanner.h"
#include "../Common/paramfile.h"
#include "../Interface/gascoigne.h"

/*---------------------------------------------------*/

namespace Gascoigne {

class P4estMeshAgentBase
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

protected:
  P4estMeshAgentBase(){};
  virtual ~P4estMeshAgentBase(){};

public:
  static std::shared_ptr<P4estMeshAgentBase> create(const ParamFile& pf);

  virtual IndexType trees_count() const = 0;
  virtual IndexType quad_count() const = 0;

  virtual void write_vtk(const std::string& fname) const = 0;
  virtual void global_refine(IndexType n = 1) = 0;
  virtual void refine_cells(IndexVector& ref) = 0;
};

} // namespace Gascoigne

/*---------------------------------------------------*/

#endif //__p4estmeshagentbase_h
