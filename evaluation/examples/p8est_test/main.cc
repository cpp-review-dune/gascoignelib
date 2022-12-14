/**
 *
 * Copyright (C) 2020 by the Gascoigne 3D authors
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

#include <Common/paramfile.h>
#include <Interface/gascoigne.h>
#include <Mesh/p4estmeshagentbase.h>

/*---------------------------------------------------*/

using namespace Gascoigne;

int
main(int argc, char** argv)
{
  // Set the name of the parameter file.
  ParamFile paramfile("config.param");

  auto pma = P4estMeshAgentBase::create(paramfile);
  IndexVector refine_cells;
  refine_cells.push_back(44);
  pma->refine_cells(refine_cells);
  pma->write_vtk("Results/solve.00000");

  return 0;
}
