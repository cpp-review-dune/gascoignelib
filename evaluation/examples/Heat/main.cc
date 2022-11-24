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
#include <Problems/functionalcontainer.h>
#include <Problems/problemcontainer.h>

#include "local.h"

#ifdef USE_CUDA
#include <cudaloop.h>
#else
#include <Solver/stdloop.h>
#endif

/*---------------------------------------------------*/

int
main(int argc, char** argv)
{
  // Set the name of the parameter file.
  Gascoigne::ParamFile paramfile("heat.param");
  if (argc >= 2) {
    paramfile.SetName(argv[1]);
  }

  Gascoigne::HeatProblem LPD;
  LPD.BasicInit(paramfile);

  Gascoigne::ProblemContainer PC;
  PC.AddProblem("heat", &LPD);

  Gascoigne::FunctionalContainer FC;

#ifdef USE_CUDA
  Gascoigne::CudaLoop loop;
#else
  Gascoigne::BasicLoop loop;
#endif
  loop.BasicInit(paramfile, &PC, &FC);
  loop.timerun("heat");

  return 0;
}
