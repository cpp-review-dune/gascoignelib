#ifndef __WeightedPointFunctional_h
#define __WeightedPointFunctional_h

#include "pointfunctional.h"

/**********************************************************/
namespace Gascoigne
{
class WeightedPointFunctional : public PointFunctional
{
  protected:
  std::vector<double> _weights;

  public:
    WeightedPointFunctional() : PointFunctional() {}
    ~WeightedPointFunctional() { }

    void BasicInit(const std::vector<Vertex2d>& v2d, const std::vector<int>& comps, const std::vector<double>& weights);
    void BasicInit(const std::vector<Vertex3d>& v3d, const std::vector<int>& comps, const std::vector<double>& weights);

    const std::vector<double>& GetWeights()    const { return _weights;}

    double J(const std::vector<double>& u) const;
  
    std::string GetName() const {return "WeightedPointFunctional";}
};
}
/**********************************************************/

#endif
