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
*/

#include <problemcontainer.h>

#include  "navierstokesproblem.h"
#include  "loop.h"
#include  "functionals.h"



/*---------------------------------------------------*/
int main(int argc, char** argv)
{
  std::string config = "bench2d1";
  
  
  /// Sets the name of the parameter file. 
  Gascoigne::ParamFile paramfile("bench2d1.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  /**
   * The problem descriptor collects all information such as equation,
   * boundary data, right hand side that are required to properly
   * describe one pde problem
   */
  Gascoigne::NavierStokesStationary<2> LPD;
  LPD.BasicInit(paramfile);
  const Gascoigne::NavierStokesData& data = LPD.GetData();

  /**
   * The ProblemContainer collects all Problems that are relevant. Here
   * we only add one problem, the LaplaceProblem specified above. Each
   * Problem is added with a label which allows us to switch between 
   * different problems.
   */
  Gascoigne::ProblemContainer PC;
  PC.AddProblem(config, &LPD);
  
  /**
   * The FunctionalContainer can collect various functionals for 
   * postprocessing of the solution, e.g. by computing point values
   * or averages of the solution. 
   */
  Gascoigne::FunctionalContainer FC;

  // squared l2-norm of the divergence
  Divergence<2> divfunctional;
  FC.AddFunctional("0 divergence",&divfunctional);
  // drag & lift
  Drag dragfunctional(config);
  Lift liftfunctional(config);
  FC.AddFunctional("1 drag",&dragfunctional);
  FC.AddFunctional("2 lift",&liftfunctional);
  // pressure difference
  PDiff pdiff(config);
  //  FC.AddFunctional("3 pdiff",&pdiff);
  // drag & Lift as boundary functional for comparison. The evaluation
  // via the Babuska-Miller-Trick is superior.
  BoundaryDrag bdragfunctional(data,config);
  BoundaryLift bliftfunctional(data,config);
  FC.AddFunctional("4 b-drag",&bdragfunctional);
  FC.AddFunctional("5 b-lift",&bliftfunctional);
  
  /**
   * The Loop controls the high level algorithmic process. Here, we
   * setup the problem, we manage the triangulation and we call the
   * different numerical method for solving the problem.
   */
  Loop loop;
  loop.BasicInit(config,paramfile, &PC, &FC);
  loop.stationaryrun(config);

  return 0;
}
