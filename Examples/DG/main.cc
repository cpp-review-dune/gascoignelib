/**
 *
 * Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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

#include "dgsolver.h"
#include "local.h"
#include "stdloop.h"

using namespace Gascoigne;

class MyLoop : public StdLoop {
public:
  void BasicInit(const ParamFile *paramfile, const ProblemContainer *PC,
                 const FunctionalContainer *FC = NULL) {
    GetMultiLevelSolverPointer() = new DGMultiLevelSolver;
    StdLoop::BasicInit(paramfile, PC, FC);
  }
};

/*---------------------------------------------------*/

int main(int argc, char **argv) {
  ParamFile paramfile("dg.param");
  if (argc >= 2) {
    paramfile.SetName(argv[1]);
  }

  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("dg", &LPD);

  FunctionalContainer FC;

  /////////////
  // Loop
  /////////////
  MyLoop loop;
  loop.BasicInit(&paramfile, &PC, &FC);

  loop.run("dg");

  return 0;
}
