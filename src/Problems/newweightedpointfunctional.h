#ifndef __NewWeightedPointFunctional_h
#define __NewWeightedPointFunctional_h

#include "newpointfunctional.h"

/**********************************************************/
namespace Gascoigne
{
class NewWeightedPointFunctional : public NewPointFunctional
{
  protected:
  std::vector<double> _weights;

  public:
    NewWeightedPointFunctional() : NewPointFunctional() {beautifulname = "NewWeightedPointFunctional"; }
    ~NewWeightedPointFunctional() { }

    void BasicInit(const std::vector<Vertex2d>& v2d, const std::vector<int>& comps, const std::vector<double>& weights);
    void BasicInit(const std::vector<Vertex3d>& v3d, const std::vector<int>& comps, const std::vector<double>& weights);

    const std::vector<double>& GetWeights()    const { return _weights;}

    double J(const std::vector<double>& u) const;
};
}
/**********************************************************/

#endif
