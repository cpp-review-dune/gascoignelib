#include "weightedpointfunctional.h"

using namespace std;

/**********************************************************/
namespace Gascoigne
{
void WeightedPointFunctional::BasicInit(const vector<Vertex2d>& v2d, const vector<int>& comps, const vector<double>& weights)
{
  _weights = weights;
  PointFunctional::BasicInit(v2d,comps);
}

/**********************************************************/

void WeightedPointFunctional::BasicInit(const vector<Vertex3d>& v3d, const vector<int>& comps, const vector<double>& weights)
{
  _weights = weights;
  PointFunctional::BasicInit(v3d,comps);
}

/**********************************************************/

double WeightedPointFunctional::J(const vector<double>& u) const
{
  double s=0;
  for(int i=0;i<u.size();++i) s+=_weights[i]*u[i];

  return s;
}
}
/**********************************************************/

