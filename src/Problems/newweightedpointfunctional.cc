#include "newweightedpointfunctional.h"

using namespace std;

/**********************************************************/
namespace Gascoigne
{
void NewWeightedPointFunctional::BasicInit(const std::vector<Vertex2d>& v2d, const std::vector<int>& comps, const std::vector<double>& weights)
{
  _weights = weights;
  NewPointFunctional::BasicInit(v2d,comps);
}

/**********************************************************/

void NewWeightedPointFunctional::BasicInit(const std::vector<Vertex3d>& v3d, const std::vector<int>& comps, const std::vector<double>& weights)
{
  _weights = weights;
  NewPointFunctional::BasicInit(v3d,comps);
}

/**********************************************************/

double NewWeightedPointFunctional::J(const std::vector<double>& u) const
{
  double s=0;
  for(int i=0;i<u.size();++i) s+=_weights[i]*u[i];

  return s;
}
}
/**********************************************************/

