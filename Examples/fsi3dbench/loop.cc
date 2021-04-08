#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"

using namespace std;

double __DT, __TIME, __THETA;

namespace Gascoigne {

template <int DIM> void Loop<DIM>::run(const std::string &problemlabel) {
  double STOP_TIME;
  int _initial_refine = 0;

  if (1) {
    DataFormatHandler DFH;
    DFH.insert("start_time", &__TIME, 0.0);
    DFH.insert("stop_time", &STOP_TIME, 0.0);
    DFH.insert("theta", &__THETA, 0.0);
    DFH.insert("dt", &__DT, 0.0);
    FileScanner FS(DFH, _paramfile, "Equation");
    assert(STOP_TIME > __TIME);
    assert(__DT > 0);
    assert(__THETA > 0);
  }
  if (1) {
    DataFormatHandler DFH;
    DFH.insert("initialrefine", &_initial_refine, 0);
    FileScanner FS(DFH, _paramfile, "Loop");
  }

  _niter = static_cast<int>((STOP_TIME - __TIME + 1.e-12) / __DT);

  VectorInterface u("u"), f("f"), old("old"), def("def"), defold("defold");

  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(old);
  GetMultiLevelSolver()->ReInitVector(f);
  GetMultiLevelSolver()->ReInitVector(def);
  GetMultiLevelSolver()->ReInitVector(defold);

  InitSolution(u);

  for (int i = 0; i < _initial_refine; ++i) {
    GlobalVector ualt;
    CopyVector(ualt, u);
    GetMeshAgent()->global_refine(1);

    GetMultiLevelSolver()->ReInit(problemlabel);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(f);
    GetMultiLevelSolver()->ReInitVector(def);
    GetMultiLevelSolver()->ReInitVector(defold);
    GetMultiLevelSolver()->InterpolateSolution(u, ualt);
  }

  GetMultiLevelSolver()->Zero(def);

  ofstream func_log("functional.txt");
  func_log.precision(12);

  // damp at 1, 2, a.s.o.
  int DI = static_cast<int>(1.0 / __DT);
  cout << DI << endl;

  StopWatch _all;

  for (_iter = 1; _iter <= _niter; _iter++) {
    cout << "========== " << _iter << ": " << __TIME << " -> " << __TIME + __DT
         << "  (" << __THETA << ")" << endl;
    __TIME += __DT;

    GetMultiLevelSolver()->Equ(old, 1.0, u);
    GetMultiLevelSolver()->Equ(defold, 1.0, def);

    GetMultiLevelSolver()->AddNodeVector("OLD", old);
    GetMultiLevelSolver()->AddNodeVector("DEFOLD", defold);
    GetMultiLevelSolver()->AddNodeVector("DEF", def);

    //	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    // assert(Solve(u,f)=="converged");
    //// SOLVE
    GetMultiLevelSolver()->GetSolver()->Zero(f);
    GetMultiLevelSolver()->GetSolver()->Rhs(f);
    GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
    // set offset first, so nodes that are both periodic and dirichlet will
    // become dirichlet
    GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
    GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
    dynamic_cast<FSIMultiLevelSolver<DIM> *>(GetMultiLevelSolver())
        ->UpdateDeformation(u);

    string solved =
        GetMultiLevelSolver()->Solve(u, f, GetSolverInfos()->GetNLInfo());
    assert(solved == "converged");
    Output(u, "Results/");

    if (_iter % 5 == 0) {
      GetMultiLevelSolver()->GetSolver()->Visu("Results/u", u, _iter);
      WriteMeshAndSolution("Results/u", u);
    }

    DoubleVector juh = Functionals(u, f);
    GetMultiLevelSolver()->DeleteNodeVector("OLD");
    GetMultiLevelSolver()->DeleteNodeVector("DEFOLD");
    GetMultiLevelSolver()->DeleteNodeVector("DEF");

    func_log << __TIME << "\t" << 0.0 << "\t" << juh << endl;
  }

  func_log.close();
}

template class Loop<2>;
template class Loop<3>;

} // namespace Gascoigne
