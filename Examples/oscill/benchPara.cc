
#include "local.h"
#include "paraLoop.h"
#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include "dirichletdatabycolor.h"

using namespace Gascoigne;
using namespace std;

bool parseCommandlineArguments(int argc, char* argv[], int& maxIterations, int& noTests,
                               std::vector<double>& coarse_theta, std::vector<double>& dtcoarse);

int main(int argc, char* argv[])
{
    std::vector<double> coarse_theta;
    std::vector<double> dtcoarse;
    int noTests, maxIterations;

    bool parse =
      parseCommandlineArguments(argc, argv, maxIterations, noTests, coarse_theta, dtcoarse);
    if (!parse)
    {
        return 1;
    }
    else
    {
        parareal<2, log_level::results>::compareParaSerial(maxIterations, dtcoarse, coarse_theta);
    }
}

bool parseCommandlineArguments(int argc, char* argv[], int& maxIterations, int& noTests,
                               std::vector<double>& coarse_theta, std::vector<double>& dtcoarse)
{
    if (argc >= 2)
    {
        int cl_mI = std::atoi(argv[1]);

        if (cl_mI > 0)
        {
            maxIterations = cl_mI;
        }
        else
        {
            return false;
        }
    }
    if (argc >= 3)
    {
        int cl_noT = std::atoi(argv[2]);

        if (cl_noT > 0)
        {
            noTests = cl_noT;
        }
        else
        {
            return false;
        }
        assert(argc >= 2 * noTests + 3);
        dtcoarse.reserve(noTests);
        coarse_theta.reserve(noTests);

        for (int i = 1; i <= noTests; i++)
        {
            dtcoarse.push_back(std::atof(argv[2 * i + 1]));
            coarse_theta.push_back(std::atof(argv[2 * (i + 1)]));
        }
    }
    return true;
}
