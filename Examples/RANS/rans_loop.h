#include "stdloop.h"
#include "problem.h"
#include "gascoignemesh2d.h"

class rans_loop : public StdLoop
{
public:
    void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,
                   const FunctionalContainer* FC)
    {
        GetMeshAgentPointer() = new MeshAgent;
        StdLoop::BasicInit(paramfile, PC, FC);
        DataFormatHandler DFH;
        DFH.insert("dt", &_dt, 0.0);
        FileScanner FS(DFH, _paramfile, "Loop");
    }

    void run(const std::string& problemlabel)
    {
        GetMultiLevelSolver()->ReInit(problemlabel);  // crates hierarchy of meshes and solvers
        set_time(_dt, _time);

        VectorInterface u("u"), f("f"), old("old");
        GetMultiLevelSolver()->ReInitVector(u);
        GetMultiLevelSolver()->ReInitVector(old);
        GetMultiLevelSolver()->ReInitVector(f);
        InitSolution(u);
        std::cout << GetMeshAgent()->ncells() << " cells" << '\n';
        GetMultiLevelSolver()->GetSolver()->OutputSettings();
        // here you can prerefine
        GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
        // TIME LOOP
        for (_iter = 1; _iter <= _niter; _iter++)
        {
            cout << "time step " << _time << " -> " << _time + _dt << endl;
            _time += _dt;
            set_time(_dt, _time);
            GetMultiLevelSolver()->Equ(old, 1.0, u);

            // Solve RANS
            GetMultiLevelSolver()->AddNodeVector("old", old);
            assert(Solve(u, f) == "converged");
            // DoubleVector juh = Functionals(u,f);
            GetMultiLevelSolver()->DeleteNodeVector("old");
        }
    }

private:
    double _dt;
    double _time;

    void set_time(const double dt, const double time) const
    {
        for (auto i = 0; i < GetMultiLevelSolver()->nlevels(); ++i)
        {
            GetMultiLevelSolver()->GetSolver(i)->GetProblemDescriptor()->SetTime(time, dt);
        }
    }
};
