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
  IndexVector refine_cells = { 1 };
  pma->refine_cells(refine_cells);

  pma->write_vtk("Results/mesh");

  auto dof = pma->create_dofhandler(1);

  GhostVectorAgent gva;
  gva.Register("u");
  gva["u"] = new GlobalVector(dof->nnodes(), 1, 0);
  (*gva["u"])[dof->nnodes() - 1] = 1.0;

  gva.Register("v");
  gva["v"] = new GlobalVector(dof->nnodes(), 1, 0);
  (*gva["v"])[dof->nnodes() - 2] = 1.0;

  gva.Register("x");
  gva["x"] = new GlobalVector(dof->nnodes(), 1, 0);
  for (auto i : dof->get_nodes_of_cell(3)) {
    (*gva["x"])[i] = 1;
  }

  gva.Register("w");
  gva["w"] = new GlobalVector(dof->nnodes(), 1, 0);
  for (IndexType i = dof->nnodes() - dof->nhanging(); i < dof->nnodes(); ++i) {
    (*gva["w"])[i] = 1;
  }

  dof->write_vtk("Results/solve.00000.vtk",
                 .0,
                 gva,
                 std::vector<std::string>({ "u", "v", "w", "x" }));

  return 0;
}
