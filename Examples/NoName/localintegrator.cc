#include  "localintegrator.h"


using namespace std;

/* ----------------------------------------- */

LocalIntegrator::LocalIntegrator() : GalerkinIntegrator<2>()
{
  cerr << "\tLocalIntegrator::LocalIntegrator()\n";
}
