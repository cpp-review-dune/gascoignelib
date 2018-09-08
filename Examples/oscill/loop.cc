#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"

using namespace std;

double __DT, __TIME, __THETA;

namespace Gascoigne
{
template <int DIM>
void Loop<DIM>::run(const std::string& problemlabel)
{
    double STOP_TIME;

    if (1)
    {
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

    _niter            = static_cast<int>((STOP_TIME - __TIME + 1.e-12) / __DT);
    auto sub_int_time = static_cast<int>(_niter / 8);
    VectorInterface u("u"), f("f"), old("old");

    GetMultiLevelSolver()->ReInit(problemlabel);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(f);
    InitSolution(u);

    ofstream func_log("functional.txt");
    func_log.precision(12);

    StopWatch _all;

    for (_iter = 1; _iter <= _niter; _iter++)
    {
        cout << "========== " << _iter << ": " << __TIME << " -> " << __TIME + __DT << "  ("
             << __THETA << ")" << endl;
        __TIME += __DT;

        for (int i = 0; i < GetMultiLevelSolver()->nlevels(); i++)
        {
            GetMultiLevelSolver()->GetSolver(i)->GetProblemDescriptor()->SetTime(__TIME, __DT);
        }

        GetMultiLevelSolver()->Equ(old, 1.0, u);
        GetMultiLevelSolver()->AddNodeVector("old", old);
        assert(Solve(u, f) == "converged");
        // GetMultiLevelSolver()->DeleteNodeVector("old");
        // GetMultiLevelSolver()->SetProblem("div");
        //
        // std::cout << "set div" << '\n';
        //
        // GetMultiLevelSolver()->Equ(old, 1.0, u);
        // GetMultiLevelSolver()->AddNodeVector("old", old);
        // assert(StdLoop::Solve(u, f) == "converged");
        // std::cout << "solved div" << '\n';
        // GetMultiLevelSolver()->SetProblem("fsi");
        if (_iter % sub_int_time == 0)
        {
            auto idx = static_cast<int>(_iter/sub_int_time);
            GetMultiLevelSolver()->GetSolver()->Visu("Results/u", u, idx);
            //WriteMeshAndSolution("Results/u", u);
            GetMultiLevelSolver()->GetSolver()->Write(u, "Results/u." + to_string(idx));
        }

        DoubleVector juh = Functionals(u, f);
        GetMultiLevelSolver()->DeleteNodeVector("old");

        func_log << __TIME << "\t" << juh << endl;
    }

    func_log.close();
}

template class Loop<2>;
template class Loop<3>;

}  // namespace Gascoigne
