#ifndef __NewWeightedDiracRightHandSide_h
#define __NewWeightedDiracRightHandSide_h

#include "newdiracrighthandside.h"
#include "newweightedpointfunctional.h"

/**********************************************************/
namespace Gascoigne
{
class NewWeightedDiracRightHandSide : public NewDiracRightHandSide
{
  protected:
  std::vector<double> _weights;

  public:
    NewWeightedDiracRightHandSide() : NewDiracRightHandSide() { }
    ~NewWeightedDiracRightHandSide() { }

    void BasicInit(const NewWeightedPointFunctional* WPF);
     
    double operator()(int i, const Vertex2d& v) const;
    double operator()(int i, const Vertex3d& v) const; 
    std::string GetName() const { return "NewWeightedDiracRightHandSide"; }
};
}
/**********************************************************/

#endif
