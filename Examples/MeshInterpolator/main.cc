#include  "projectiondescriptor.h"
#include  "stdmultilevelsolver.h"
#include  "meshinterpolator.h"
#include  "solverinfos.h"
#include  "meshagent.h"
#include  "backup.h"

using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  string finename = "start2";
  
  ProjectionProblemDescriptor PD;
  PD.BasicInit(&paramfile);
  
  ///////////////////////
  // MESH 
  ///////////////////////

  MeshAgent MA;
  MA.BasicInit(&paramfile);

  ///////////////////////
  // MG Solver
  ///////////////////////

  StdMultiLevelSolver MLS;
  MLS.BasicInit(&MA,&paramfile);
  MLS.ReInit(PD);
  MLS.GetSolver()->OutputSettings();
    
  SolverInfos SI;
  SI.BasicInit(&paramfile);
  SI.GetNLInfo().control().matrixmustbebuild() = 1;

  ///////////////////////
  // Vectors
  ///////////////////////

  VectorInterface u("u"), f("f");
  MLS.ReInitVector(u);
  MLS.ReInitVector(f);
  MLS.ReInit(PD);
  
  MLS.Zero(u);
  MLS.Zero(f);
 
  ///////////////////////
  // Mesh Interpolator
  ///////////////////////
    {
      MeshInterpolator MI;
      MI.BasicInit(MLS.GetSolver()->GetDiscretization(),&MA,finename);
      
      GlobalVector uold;
      ReadBackUpResize(uold,finename+".bup");
      MI.RhsForProjection(MLS.GetSolver()->GetGV(f),uold);
    }
  string coarsename = "Results/" + finename + ".projected";

  Monitor moni;
  moni.init(&paramfile,1);
  moni.set_directory("Results");
  MLS.SetMonitorPtr(&moni);

  MLS.Solve(MLS.nlevels()-1,u,f,SI.GetNLInfo());
  MLS.GetSolver()->Visu(coarsename,u,0);

  return 0;
}
