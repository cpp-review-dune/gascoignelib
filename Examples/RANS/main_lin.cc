#include "problem.h"
#include <meshagent.h>
#include "problemdescriptorbase.h"
#include "rans_loop.h"

using namespace Gascoigne;

int main(int argc, char** argv)
{
    ParamFile paramfile("rans.param");

    LinProblemDescriptor rans_pd;
    rans_pd.BasicInit(&paramfile);

    ProblemContainer PC;
    PC.AddProblem("rans", &rans_pd);

    rans_loop loop;
    loop.BasicInit(&paramfile, &PC, nullptr);
    loop.run("rans");

    return 0;
}
