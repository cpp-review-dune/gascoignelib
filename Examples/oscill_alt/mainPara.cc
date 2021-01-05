#include "local.h"
#include "paraLoop.h"
#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include "dirichletdatabycolor.h"

using namespace Gascoigne;
using namespace std;

//--------------------------------------------------------------------------------------------------
bool parseCommandlineArguments(int     argc,
                               char*   argv[],
                               int&    maxIterations,
                               double& coarse_theta,
                               double& dtcoarse,
                               double& fine_theta,
                               double& dtfine);
//--------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
#ifdef BDF
  parareal<2>::paramfile= ParamFile("fsi-3.param");
#else
  parareal<2>::paramfile= ParamFile("fsi-wm.param");
#endif
  int    maxIterations;
  double coarse_theta;
  double dtcoarse;
  double fine_theta;
  double dtfine;

  bool parse= parseCommandlineArguments(
    argc, argv, maxIterations, coarse_theta, dtcoarse, fine_theta, dtfine);
  if (!parse) {
    return 1;
  }
  if (parse && argc == 6) {
    if (fine_theta == 0 && dtfine == 0) {
      parareal<2>::runPara(maxIterations, coarse_theta, dtcoarse);
    } else {
      parareal<2>::runPara(maxIterations, coarse_theta, dtcoarse, fine_theta, dtfine);
    }
  } else if (parse && argc == 2) {
    parareal<2>::runPara(maxIterations);
  } else {
    parareal<2>::runPara();
  }
}

bool parseCommandlineArguments(int     argc,
                               char*   argv[],
                               int&    maxIterations,
                               double& coarse_theta,
                               double& dtcoarse,
                               double& fine_theta,
                               double& dtfine) {
  if (argc >= 2) {
    int cl_mI= std::atoi(argv[1]);

    if (cl_mI > 0) {
      maxIterations= cl_mI;
    } else {
      return false;
    }
  }
  if (argc >= 3)  // suppose that when this is the case all arguments are given via cl
  {
    double cl_ctt= std::atof(argv[2]);
    double cl_dtc= std::atof(argv[3]);
    double cl_tt = std::atof(argv[4]);
    double cl_dtf= std::atof(argv[5]);
    if (cl_ctt > 0 && cl_dtc > 0 && cl_tt > 0 && cl_dtf > 0) {
      coarse_theta= cl_ctt;
      dtcoarse    = cl_dtc;
      fine_theta  = cl_tt;
      dtfine      = cl_dtf;
    } else {
      return false;
    }
  }

  return true;
}
