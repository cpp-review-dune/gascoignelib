#include "weighteddiracrighthandside.h"

using namespace std;

/**********************************************************/
namespace Gascoigne
{

void WeightedDiracRightHandSide::BasicInit(const WeightedPointFunctional* WPF)
{
  _v2d = WPF->GetPoints2d();
  _v3d = WPF->GetPoints3d();

  _comps = WPF->GetComps();

  _weights = WPF->GetWeights();
}

/**********************************************************/

double WeightedDiracRightHandSide::operator()(int i,const Vertex2d& v) const
{
  return _weights[i];
}

/**********************************************************/

double WeightedDiracRightHandSide::operator()(int i,const Vertex3d& v) const
{
  return _weights[i];
}

}
