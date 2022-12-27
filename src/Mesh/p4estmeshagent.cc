/**
 *
 * Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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

#include "p4estmeshagent.h"

#include "p4estmeshagent2d.h"
#include "p4estmeshagent3d.h"

namespace Gascoigne {

std::shared_ptr<P4estMeshAgent>
P4estMeshAgent::create(const ParamFile& pf)
{
  IndexType dimension;
  IndexType prerefine;
  std::string gridname;

  DataFormatHandler DFH;
  DFH.insert("dimension", &dimension, 0);
  DFH.insert("gridname", &gridname, "none");
  DFH.insert("prerefine", &prerefine, 0);

  FileScanner FS(DFH);
  FS.readfile(pf, "Mesh");

  switch (dimension) {
    case 2:
      return std::make_shared<P4estMeshAgent2d>(gridname, prerefine);
    case 3:
      return std::make_shared<P4estMeshAgent3d>(gridname, prerefine);
  }
  return nullptr;
}
} // namespace Gascoigne