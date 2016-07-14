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


#include  "local.h"
#include "stopwatch.h"

extern Gascoigne::Stoppers GlobalStopWatch;

using namespace std;
using namespace Gascoigne;

double __DT;
double __TIME;


/*---------------------------------------------------*/

class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC) 
    {
      GetMeshAgentPointer() = new BenchMarkMeshAgent;
      StdLoop::BasicInit(paramfile, PC, FC);
    }

  void run(const std::string& problemlabel)
  {
    VectorInterface u("u"), f("f"), old("old");


    GetMultiLevelSolver()->ReInit(problemlabel);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(f);
    GetMultiLevelSolver()->GetSolver()->OutputSettings();
    InitSolution(u);
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

    PrintMeshInformation();

    double stoptime;
    DataFormatHandler DFH;
    DFH.insert("dt" , &__DT , 0.);
    DFH.insert("time",&__TIME, 0.);
    DFH.insert("stoptime",&stoptime, 0.);
    FileScanner FS(DFH, _paramfile, "Loop");
    assert(__DT>0);
    
    
    for (_iter=1; __TIME+1.e-10 < stoptime; _iter++)
      {
	cout << "\n  ============== " 
	     << __TIME << " -> " << __TIME + __DT << endl;
	__TIME+=__DT;
	
	GetMultiLevelSolver()->Equ(old,1.0,u);

	GetMultiLevelSolver()->AddNodeVector("old",old);
	Solve(u,f);
	DoubleVector juh = Functionals(u,f);
	GetMultiLevelSolver()->DeleteNodeVector("old");


	
	_clock_functionals.stop();
      }
  }
  
  
  
};

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  GlobalStopWatch.start("0 Alles\t");
  
  ParamFile paramfile("bench.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);
  
  ProblemContainer PC;
  PC.AddProblem("navier stokes", &LPD);

	    // Functionals
  FunctionalContainer FC;
  LocalDomainFunctionals_FlagForce df_flagforce_lift("lift");
  LocalDomainFunctionals_FlagForce df_flagforce_drag("drag");
  FC.AddFunctional("lift",&df_flagforce_lift);
  FC.AddFunctional("drag",&df_flagforce_drag);
  
  LocalLoop loop;
  loop.BasicInit(&paramfile, &PC, &FC);

  loop.run("navier stokes");

  GlobalStopWatch.stop("0 Alles\t");

  GlobalStopWatch.print();

  
  
  return 0;
}