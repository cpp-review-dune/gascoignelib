#ifndef __WeightedDiracRightHandSide_h
#define __WeightedDiracRightHandSide_h

#include "diracrighthandside.h"
#include "weightedpointfunctional.h"

/**********************************************************/
namespace Gascoigne
{
class WeightedDiracRightHandSide : public DiracRightHandSide
{
  protected:
  std::vector<double> _weights;

  public:
    WeightedDiracRightHandSide() : DiracRightHandSide() { }
    ~WeightedDiracRightHandSide() { }

    void BasicInit(const WeightedPointFunctional* WPF);
     
    double operator()(int i, const Vertex2d& v) const;
    double operator()(int i, const Vertex3d& v) const; 
    std::string GetName() const { return "WeightedDiracRightHandSide"; }
};
}
/**********************************************************/

#endif
