#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "localhierarchicalmesh3d.h"
#include "localmeshagent.h"
#include <time.h>
using namespace std;

namespace Gascoigne {

template<int DIM>
void
Loop<DIM>::run(const std::string& problemlabel)
{
  ofstream func_log("functional.txt");
  func_log.precision(12);
  func_log << "cells\t"
           << "nodes\t"
           << "func\t\t\t"
           << "newton\t"
           << "GMRES\t"
           << "GMRES_mean\t"
           << "time" << endl;
  clock_t t;

  Vector u("u"), f("f");

  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(f);

  int refinem_steps = 1;
  for (int i = 1; i <= refinem_steps; i++) {

    StopWatch _all;

    t = clock();
    cout << "=================================================" << endl;
    cout << "Computation on preref:" << i << endl;
    //// SOLVE

    GetMultiLevelSolver()->GetSolver()->Zero(u);
    GetMultiLevelSolver()->GetSolver()->Zero(f);
    GetMultiLevelSolver()->GetSolver()->Rhs(f);
    GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
    // set offset first, so nodes that are both periodic and dirichlet will
    // become dirichlet
    GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
    GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
    GetSolverInfos()->GetNLInfo().statistics().reset();
    GetSolverInfos()->GetNLInfo().GetLinearInfo().statistics().reset();

    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    string solved =
      GetMultiLevelSolver()->Solve(u, f, GetSolverInfos()->GetNLInfo());

    NLInfo& nlinfo = GetSolverInfos()->GetNLInfo();
    // anzahl newtonschritte
    double newton_Steps = nlinfo.statistics().totaliter();
    cout << "newton_Steps: " << newton_Steps << endl;
    // GMRES Iter(ausgabe gibt bei 0 Iter [1]!!)
    double GMRES_Steps = nlinfo.GetLinearInfo().statistics().totaliter();
    cout << "GMRES_Steps: " << GMRES_Steps << endl;
    cout << "GMRES_mean: " << GMRES_Steps / newton_Steps << endl;
    // anzahl matrix aufbauen(erste nicht mitgezaehlt!)"
    cout << "matrix build: " << nlinfo.statistics().totalmatrix() << endl;

    assert(solved == "converged");
    DoubleVector juh = Functionals(u, f);

    t = clock() - t;

    func_log << GetMeshAgent()->ncells() << "\t" << GetMeshAgent()->nnodes()
             << "\t" << juh << "\t" << newton_Steps << "\t" << GMRES_Steps
             << "\t" << GMRES_Steps / newton_Steps << "\t\t"
             << ((float)t) / CLOCKS_PER_SEC << endl;
    GetMultiLevelSolver()->GetSolver()->Visu("Results/u", u, i);

    // interpolate solution & refine
    if (i < refinem_steps) {
      GetMeshAgent()->global_refine(1);
      GetMultiLevelSolver()->ReInit(problemlabel);
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(f);
      GetMultiLevelSolver()->Zero(u);
    }
  }
  func_log.close();

  string name = "Results/transformed_mesh";
  LocalMeshAgent* lMA = dynamic_cast<LocalMeshAgent*>(GetMeshAgent());
  lMA->write_gup(name, GetMultiLevelSolver()->GetSolver()->GetGV(u));
  cout << " " << name << ".gup" << endl;

  GetMultiLevelSolver()->DeleteNodeVector("u");
  GetMultiLevelSolver()->DeleteNodeVector("f");
}

template class Loop<2>;
template class Loop<3>;

} // namespace Gascoigne
