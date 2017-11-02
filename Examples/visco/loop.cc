#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include  "compose_name.h"
using namespace std;


double __DT,__TIME;


//extern bool BOUNDARY;


namespace Gascoigne
{

  template<int DIM>
  void Loop<DIM>::run(const std::string& problemlabel)
  {
    double STOP_TIME;
    int _initial_refine=0;

    if (1)
      {
  DataFormatHandler DFH;
  DFH.insert("start_time" ,    &__TIME , 0.0);
  DFH.insert("stop_time" ,    &STOP_TIME , 0.0);
  DFH.insert("dt" ,    &__DT,0.0);
  FileScanner FS(DFH, _paramfile, "Equation");
  assert(STOP_TIME>__TIME);
  assert(__DT>0);
      }
    if (1)
      {
  DataFormatHandler DFH;
  DFH.insert("initialrefine", &_initial_refine,0);
  FileScanner FS(DFH, _paramfile, "Loop");
      }



    _niter = static_cast<int> ( (STOP_TIME-__TIME+1.e-12)/__DT );

    VectorInterface u("u"), f("f"),  old("old");

    GetMultiLevelSolver()->ReInit(problemlabel);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(f);

    GetMultiLevelSolver()->SetProblem("vps");
    InitSolution(u);

    ofstream func_log("functional.txt");
    func_log.precision(12);

    GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,0);

    double writeevery=0.01;
    double writeiter =0.0;

    for (_iter=1; _iter<=_niter; _iter++)
    {

      cout << "================================ " << _iter << ": " << __TIME << " -> "
           << __TIME+__DT <<  endl;
      __TIME += __DT;
      writeiter+=__DT;

      GetMultiLevelSolver()->Equ(old,1.0,u);

      //VPS step
      GetMultiLevelSolver()->SetProblem("vps"); BOUNDARY = true;      //this name referest to what?
//      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

      GetMultiLevelSolver()->AddNodeVector("old",old);
      assert(Solve(u,f)=="converged");
      GetMultiLevelSolver()->DeleteNodeVector("old");



      if (writeiter+1.e-8>writeevery){

        GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,_iter);

        string name;
        name = "Results/u";compose_name(name,_iter);
        GetMultiLevelSolver()->GetSolver()->Write(u,name);

        writeiter=0;
      }

      nvector<double> jkin =   Functionals(u,f);
      nvector<double> jel =   Functionals(u,f);

      func_log << __TIME << "\t" << jkin[0] << "\t" << jel[1] << endl;
    }

    func_log.close();
  }


  template class Loop<2>;
  template class Loop<3>;


}

