#include "newweighteddiracrighthandside.h"

using namespace std;

/**********************************************************/
namespace Gascoigne
{

void NewWeightedDiracRightHandSide::BasicInit(const NewWeightedPointFunctional* WPF)
{
  _v2d = WPF->GetPoints2d();
  _v3d = WPF->GetPoints3d();

  _comps = WPF->GetComps();

  _weights = WPF->GetWeights();
}

/**********************************************************/

double NewWeightedDiracRightHandSide::operator()(int i,const Vertex2d& v) const
{
  return _weights[i];
}

/**********************************************************/

double NewWeightedDiracRightHandSide::operator()(int i,const Vertex3d& v) const
{
  return _weights[i];
}

}
