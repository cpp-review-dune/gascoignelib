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

#include <memory>
#include <vector>

#include <Common/paramfile.h>
#include <Interface/gascoigne.h>
#include <Mesh/p4estmeshagent.h>
#include <Solver/ghostagent.h>

/*---------------------------------------------------*/

using namespace Gascoigne;

int
main(int argc, char** argv)
{
  // Set the name of the parameter file.
  ParamFile paramfile("config.param");

  auto pma = P4estMeshAgent::create(paramfile);
  IndexVector refine_cells = { 44 };
  pma->refine_cells(refine_cells);

  auto dof = pma->create_dofhandler(2);

  GhostVectorAgent gva;
  gva.Register("u");
  gva["u"] = new GlobalVector(dof->num_nodes(), 1, 0);
  (*gva["u"])[13] = 1.0;

  gva.Register("v");
  gva["v"] = new GlobalVector(dof->num_nodes(), 1, 0);
  (*gva["v"])[19] = 1.0;

  dof->write_vtk(
    "Results/solve.00000.vtk", .0, gva, std::vector<std::string>({ "u", "v" }));

  return 0;
}
