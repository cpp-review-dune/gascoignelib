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

#include <functionalcontainer.h>
#include <problemcontainer.h>

#include "local.h"
#include "paramfile.h"

#ifdef USE_CUDA
#include <cudaloop.h>
#else
#include <stdloop.h>
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

    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
  Gascoigne::ProblemContainer PC;
  PC.AddProblem("heat", &LPD);

    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
  // Functionals
  Gascoigne::FunctionalContainer FC;

    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \

#ifdef USE_CUDA
  Gascoigne::CudaLoop loop;
#else
  Gascoigne::BasicLoop loop;
#endif
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
  loop.BasicInit(paramfile, &PC, &FC);
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
  loop.timerun("heat");
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \

  return 0;
}
