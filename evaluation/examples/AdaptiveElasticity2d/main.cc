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
#include <DofHandler/dofhandler.h>
#include <Problems/functionalcontainer.h>
#include <Problems/problemcontainer.h>

#ifdef USE_CUDA
#include <cudaloop.h>
#else
#include <Solver/stdloop.h>
#endif

#include "local.h"

#define HEAT_COMP 2

double __TIME;

namespace Gascoigne {
extern Timer GlobalTimer;

#ifdef USE_CUDA
class MyLoop : public CudaLoop
#else
class MyLoop : public StdLoop
#endif
{
public:
  void adaptiverun(const std::string& problemlabel)
  {
    Vector u("u");
    Vector f("f");
    Vector old("old");
    Matrix A("A");

    GetMultiLevelSolver()->ReInit();

    // initial refinement of the mesh
    for (size_t adarefine = 0; adarefine < 2; ++adarefine) {
      const DofHandler2d* M2d = dynamic_cast<const DofHandler2d*>(
        GetMultiLevelSolver()->GetSolver()->GetMesh());
      assert(M2d);
      IndexVector ref, coa;
      for (IndexType i = 0; i < M2d->nnodes(); ++i) {
        const Vertex2d& v = M2d->vertex2d(i);
        if (((v.x() == 0) && (v.y() == 0)) || (v.x() == 2)) {
          ref.push_back(i);
        }
      }
      GetMeshAgent()->refine_nodes(ref, coa);
      GetMultiLevelSolver()->ReInit();
    }

    GetMultiLevelSolver()->SetProblem(problemlabel);
    GetMultiLevelSolver()->ReInitMatrix(A);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(f);
    GetMultiLevelSolver()->ReInitVector(old);
    InitSolution(u);

    double dt;
    DataFormatHandler DFH;
    DFH.insert("dt", &dt, 1.);
    FileScanner FS(DFH, _paramfile, "Equation");
    assert(dt > 0);

    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

    PrintMeshInformation();
    GlobalTimer.reset();

    __TIME = 0;
    for (_iter = 1; _iter <= _niter; _iter++) {
      GlobalTimer.start("iteration");

      std::cout << "\n******************** Time Step " << _iter << "\t"
                << __TIME << " -> " << __TIME + dt << std::endl;
      __TIME += dt;
      GetMultiLevelSolver()->Equ(old, 1., u);

      GetMultiLevelSolver()->AddNodeVector("OLD", old);
      Solve(A, u, f);
      GetMultiLevelSolver()->DeleteNodeVector("OLD");
      GlobalTimer.stop("iteration");
      // GlobalTimer.print();
    }
    GlobalTimer.print();
  }
};

}

/*---------------------------------------------------*/

int
main(int argc, char** argv)
{

  // Set the name of the parameter file.
  Gascoigne::ParamFile paramfile("config.param");
  if (argc >= 2) {
    paramfile.SetName(argv[1]);
  }

  // <**>
  Gascoigne::HeatProblem<HEAT_COMP> LPD;
  // <**>

  LPD.BasicInit(paramfile);

  Gascoigne::ProblemContainer PC;
  PC.AddProblem("heat", &LPD);

  // Functionals
  Gascoigne::FunctionalContainer FC;
  Gascoigne::MyLoop loop;
  loop.BasicInit(paramfile, &PC, &FC);
  loop.adaptiverun("heat");

  return 0;
}

template class Gascoigne::HeatRhs<HEAT_COMP>;
