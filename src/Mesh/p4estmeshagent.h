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

#ifndef __p4estmeshagent_h
#define __p4estmeshagent_h

#include <memory>
#include <string>

#include "../Common/dataformathandler.h"
#include "../Common/filescanner.h"
#include "../Common/paramfile.h"
#include "../DofHandler/p4estdofhandler.h"
#include "../Interface/gascoigne.h"

/*---------------------------------------------------*/

namespace Gascoigne {

class P4estMeshAgent
{
public:
  struct pquadrant_data_t
  {
    bool refine; ///< flag when true gets refinde in refine_cells
  };

protected:
  P4estMeshAgent(){};
  virtual ~P4estMeshAgent(){};

public:
  static std::shared_ptr<P4estMeshAgent> create(const ParamFile& pf);

  virtual IndexType num_cells() const = 0;

  /** Write the mesh as vtk without data */
  virtual void write_vtk(const std::string& fname) const = 0;
  virtual void global_refine(IndexType n = 1) = 0;
  virtual void refine_cells(IndexVector& ref) = 0;

  virtual std::shared_ptr<P4estDofHandler> create_dofhandler(
    IndexType degree) const = 0;
};

} // namespace Gascoigne

/*---------------------------------------------------*/

#endif //__p4estmeshagent_h
